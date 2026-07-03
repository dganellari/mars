"""Tiny arithmetic-expression AST + parser + C renderer.

This is the whole "language" the operator front-end needs for a pointwise
weak-form action: numbers, identifiers (deriv, dt1, dt2, g0, g1, g2, ...),
the four arithmetic ops, unary minus, and parentheses. Kept deliberately
minimal -- the durable part of the compiler is this grammar + the IR, not the
emitter. No external deps.
"""

# Plain classes (no dataclasses / no PEP-585 annotations) so the compiler runs
# on stock cluster Python 3.6+ as well as modern Python.


# --- AST ---------------------------------------------------------------

class Expr:
    def free_vars(self):
        raise NotImplementedError

    def render(self, names):
        raise NotImplementedError


class Num(Expr):
    def __init__(self, value):
        self.value = value  # text so 0.5 stays 0.5, not 0.50000001

    def free_vars(self):
        return set()

    def render(self, names):
        return self.value

    def diff(self, var):
        return Num("0")


class Var(Expr):
    def __init__(self, name):
        self.name = name

    def free_vars(self):
        return {self.name}

    def render(self, names):
        # Unknown identifiers pass through unchanged (constants defined in the
        # emitted scope, e.g. a material coefficient).
        return names.get(self.name, self.name)

    def diff(self, var):
        return Num("1") if self.name == var else Num("0")


class Unary(Expr):
    def __init__(self, op, operand):
        self.op = op            # '-'
        self.operand = operand

    def free_vars(self):
        return self.operand.free_vars()

    def render(self, names):
        return f"(-{self.operand.render(names)})"

    def diff(self, var):
        return Unary("-", self.operand.diff(var))


class Binary(Expr):
    def __init__(self, op, lhs, rhs):
        self.op = op            # + - * /
        self.lhs = lhs
        self.rhs = rhs

    def free_vars(self):
        return self.lhs.free_vars() | self.rhs.free_vars()

    def render(self, names):
        return f"({self.lhs.render(names)} {self.op} {self.rhs.render(names)})"

    def diff(self, var):
        a, b = self.lhs, self.rhs
        da, db = a.diff(var), b.diff(var)
        if self.op == "+":
            return Binary("+", da, db)
        if self.op == "-":
            return Binary("-", da, db)
        if self.op == "*":                       # product rule
            return Binary("+", Binary("*", da, b), Binary("*", a, db))
        # self.op == "/" : quotient rule (da*b - a*db) / (b*b)
        return Binary("/", Binary("-", Binary("*", da, b), Binary("*", a, db)),
                      Binary("*", b, b))


# --- simplification (keeps derivative expressions readable / small) ----

def _is_zero(e: Expr) -> bool:
    return isinstance(e, Num) and float(e.value) == 0.0


def _is_one(e: Expr) -> bool:
    return isinstance(e, Num) and float(e.value) == 1.0


def simplify(e: Expr) -> Expr:
    """Bottom-up identity simplification: fold 0/1 so d/dx output stays small.

    Does NOT collect like terms (g2*deriv + g2*deriv stays as-is) -- that is a
    readability nicety, not needed for correctness, and the numerical gates
    validate the result regardless.
    """
    if isinstance(e, (Num, Var)):
        return e
    if isinstance(e, Unary):
        o = simplify(e.operand)
        if _is_zero(o):
            return Num("0")
        if isinstance(o, Unary) and o.op == "-":   # -(-x) = x
            return o.operand
        return Unary("-", o)
    a, b, op = simplify(e.lhs), simplify(e.rhs), e.op
    if op == "+":
        if _is_zero(a):
            return b
        if _is_zero(b):
            return a
    elif op == "-":
        if _is_zero(b):
            return a
        if _is_zero(a):
            return simplify(Unary("-", b))
    elif op == "*":
        if _is_zero(a) or _is_zero(b):
            return Num("0")
        if _is_one(a):
            return b
        if _is_one(b):
            return a
    elif op == "/":
        if _is_zero(a):
            return Num("0")
        if _is_one(b):
            return a
    return Binary(op, a, b)


# --- tokenizer ---------------------------------------------------------

def _tokenize(s):
    toks = []
    i, n = 0, len(s)
    while i < n:
        c = s[i]
        if c.isspace():
            i += 1
            continue
        if c in "+-*/()":
            toks.append(("op", c))
            i += 1
            continue
        if c.isdigit() or c == ".":
            j = i
            while j < n and (s[j].isdigit() or s[j] in ".eE" or
                             (s[j] in "+-" and s[j - 1] in "eE")):
                j += 1
            toks.append(("num", s[i:j]))
            i = j
            continue
        if c.isalpha() or c == "_":
            j = i
            while j < n and (s[j].isalnum() or s[j] == "_"):
                j += 1
            toks.append(("id", s[i:j]))
            i = j
            continue
        raise ValueError(f"unexpected character {c!r} in expression {s!r}")
    toks.append(("end", ""))
    return toks


# --- recursive-descent parser -----------------------------------------

class _Parser:
    def __init__(self, toks):
        self.toks = toks
        self.pos = 0

    def peek(self):
        return self.toks[self.pos]

    def eat(self):
        t = self.toks[self.pos]
        self.pos += 1
        return t

    def parse(self) -> Expr:
        e = self._expr()
        if self.peek()[0] != "end":
            raise ValueError(f"trailing tokens after expression: {self.peek()}")
        return e

    def _expr(self):
        node = self._term()
        while self.peek() == ("op", "+") or self.peek() == ("op", "-"):
            op = self.eat()[1]
            node = Binary(op, node, self._term())
        return node

    def _term(self):
        node = self._factor()
        while self.peek() == ("op", "*") or self.peek() == ("op", "/"):
            op = self.eat()[1]
            node = Binary(op, node, self._factor())
        return node

    def _factor(self):
        if self.peek() == ("op", "-"):
            self.eat()
            return Unary("-", self._factor())
        return self._primary()

    def _primary(self):
        kind, val = self.peek()
        if kind == "num":
            self.eat()
            return Num(val)
        if kind == "id":
            self.eat()
            return Var(val)
        if self.peek() == ("op", "("):
            self.eat()
            e = self._expr()
            if self.peek() != ("op", ")"):
                raise ValueError("missing closing paren")
            self.eat()
            return e
        raise ValueError(f"unexpected token {self.peek()}")


def parse_expr(s: str) -> Expr:
    return _Parser(_tokenize(s)).parse()
