"""Build the opt-in metadata tracer and export fixed protocol checks."""

from collections import Counter
import json
from pathlib import Path
import subprocess


def exchange_header(include, directory):
    # Use the installed declaration without pulling CUDA internals into a host tracer.
    text = (include / '_hypre_utilities.h').read_text()
    marker = '#ifndef hypre_EXCHANGE_DATA_HEADER'
    if text.count(marker) != 1:
        raise ValueError('unsupported exchange header')
    fragment, end, _ = text.split(marker, 1)[1].partition('#endif')
    if not end or '#if' in fragment or 'hypre_DataExchangeResponse' not in fragment or 'hypre_DataExchangeList(' not in fragment:
        raise ValueError('unsupported exchange declaration')
    target = directory / 'setup-protocol-exchange.h'
    target.write_text(marker + fragment + '#endif\n')
    return target


def prepare(directory, identity, build_cache, environment, public):
    from simple_hypre_defaults import cached_toolchain, find_compiler, library_kind, matching_install
    from simple_snapshot_compare import digest, require
    import simple_startup_probe as startup

    public['failed_check'] = 'protocol_existing_preload'
    require(not environment.get('LD_PRELOAD'), 'existing_preload')
    public['failed_check'] = 'protocol_headers'
    library, include = matching_install(identity['libraries'])
    public['failed_check'] = 'protocol_build_configuration'
    compiler, flags, links, inputs = cached_toolchain(build_cache, True)
    public['failed_check'] = 'protocol_compiler_unavailable'
    compiler = find_compiler(compiler)
    require(compiler is not None, 'compiler_unavailable')
    source = Path(__file__).resolve().parent.parent / 'tests/reference/openaccel/simple_performance/setup_protocol_trace.cpp'
    target = directory / 'setup-protocol.so'
    inputs.update({str(source): digest(source), str(library): digest(library)})
    inputs.update({str(p): digest(p) for p in include.glob('*.h')})
    public['failed_check'] = 'protocol_exchange_header'
    header = exchange_header(include, directory)
    inputs[str(header)] = digest(header)
    command = [compiler, '-std=c++11', '-O2', '-fPIC', '-shared', '-I' + str(include), '-I' + str(directory)]
    command += flags + [str(source)] + links + ['-ldl', '-o', str(target)]
    public['failed_check'] = 'protocol_compile'
    with (directory / 'protocol-build.log').open('xb') as log:
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, timeout=120)
    require(result.returncode == 0 and target.is_file(), 'protocol_compile')
    public['failed_check'] = 'protocol_mpi_library'
    libraries = startup.runtime_libraries(target)
    expected = {v for p, v in identity['libraries'].items() if library_kind(p) == 'mpi'}
    actual = {v for p, v in libraries.items() if library_kind(p) == 'mpi'}
    require(actual and actual <= expected, 'protocol_mpi_library')
    public['failed_check'] = 'protocol_build_inputs_changed'
    require(all(digest(p) == v for p, v in inputs.items()), 'protocol_build_inputs_changed')
    startup.write_json(directory / 'protocol-build.json', dict(command=command, inputs=inputs,
                       libraries=libraries, binary=str(target), sha256=digest(target)))
    return dict(preload=str(target), hypre=str(library))


def recheck(directory):
    from simple_snapshot_compare import digest, require
    record = json.loads((directory / 'protocol-build.json').read_text())
    require(all(digest(p) == v for p, v in dict(record['inputs'], **record['libraries']).items())
            and digest(record['binary']) == record['sha256'], 'protocol_build_inputs_changed')


def integer(value):
    return type(value) is int and 0 <= value < 2**63


def read_events(path):
    if path.stat().st_size > 64 * 1024 * 1024:
        raise ValueError('trace too large')
    text = path.read_text()
    if not text.endswith('\n'):
        raise ValueError('partial final event')
    shapes = dict(setup_entry=set(), setup_return=set(), exchange_return=set(),
                  exchange_entry={'round', 'rank', 'ranks', 'comm_world'},
                  call_return={'id', 'ok'}, tested={'requests', 'ok', 'complete'})
    shapes.update({name: {'id', 'peer', 'tag', 'zero', 'same_comm'}
                   for name in ('send', 'recv', 'isend', 'irecv')})
    shapes.update({name: {'requests'} for name in ('test', 'testall', 'wait', 'waitall')})
    result = []
    for line in text.splitlines():
        value = json.loads(line)
        if (not isinstance(value, dict) or not isinstance(value.get('event'), str)
                or value['event'] not in shapes):
            raise ValueError('unsupported trace event')
        if set(value) != shapes[value['event']] | {'event', 'exchange'} or not integer(value['exchange']):
            raise ValueError('invalid trace event')
        for key in ('id', 'peer', 'tag', 'round', 'rank', 'ranks'):
            if key in value and not integer(value[key]):
                raise ValueError('invalid trace integer')
        for key in ('zero', 'same_comm', 'ok', 'complete', 'comm_world'):
            if key in value and type(value[key]) is not bool:
                raise ValueError('invalid trace flag')
        if 'requests' in value and (not isinstance(value['requests'], list)
                or len(value['requests']) > 100000 or not all(integer(v) and v > 0 for v in value['requests'])):
            raise ValueError('untracked request')
        result.append(value)
    return result


def trace_state(events, rank, ranks):
    if not events:
        raise ValueError('empty trace')
    current, active, setup = None, False, False
    for value in events:
        kind = value['event']
        if kind == 'setup_entry':
            if setup or active:
                raise ValueError('nested setup')
            setup = True
        elif kind == 'setup_return':
            if not setup or active:
                raise ValueError('invalid setup return')
            setup = False
        elif kind == 'exchange_entry':
            if (not setup or active or value['rank'] != rank or value['ranks'] != ranks
                    or value['round'] not in (1, 2)):
                raise ValueError('unsupported exchange')
            if current and value['exchange'] != current['entry']['exchange'] + 1:
                raise ValueError('missing exchange')
            current = dict(entry=value, messages={}, returned={}, wait=None, errors=False)
            active = True
        else:
            if not active or value['exchange'] != current['entry']['exchange']:
                raise ValueError('event outside active exchange')
            if kind == 'exchange_return':
                active = False
            elif kind in ('send', 'recv', 'isend', 'irecv'):
                if value['id'] <= 0 or value['id'] in current['messages'] or value['peer'] >= ranks:
                    raise ValueError('invalid message')
                current['messages'][value['id']] = value
            elif kind == 'call_return':
                if value['id'] not in current['messages'] or value['id'] in current['returned']:
                    raise ValueError('invalid call return')
                current['returned'][value['id']] = value['ok']
                current['errors'] |= not value['ok']
            else:
                if any(i not in current['messages'] or current['messages'][i]['event'] not in ('isend', 'irecv')
                       for i in value['requests']):
                    raise ValueError('request is not a posted operation')
                if kind == 'tested':
                    if current['wait'] is None or value['requests'] != current['wait']['requests']:
                        raise ValueError('completion without wait')
                    current['errors'] |= not value['ok']
                    if value['complete']:
                        current['wait'] = None
                else:
                    current['wait'] = value
    return current, active, setup


def active_exchange(events, rank, ranks):
    current, active, setup = trace_state(events, rank, ranks)
    if not setup or not active or current is None:
        raise ValueError('no unfinished exchange')
    return current


TRACE_REJECTIONS = {
    'trace too large': 'size_limit', 'partial final event': 'partial_record',
    'unsupported trace event': 'unsupported_event', 'invalid trace event': 'invalid_event',
    'invalid trace integer': 'invalid_integer', 'invalid trace flag': 'invalid_flag',
    'untracked request': 'untracked_request', 'empty trace': 'empty_trace',
    'nested setup': 'nested_setup', 'invalid setup return': 'invalid_setup_return',
    'unsupported exchange': 'unsupported_exchange', 'missing exchange': 'missing_exchange',
    'event outside active exchange': 'event_outside_exchange', 'invalid message': 'invalid_message',
    'invalid call return': 'invalid_call_return', 'request is not a posted operation': 'unposted_request',
    'completion without wait': 'completion_without_wait',
}


def rank_trace(directory, rank, ranks):
    summary = dict(records_valid=False, state='unknown', rejection='unreadable_record')
    try:
        events = read_events(directory / ('protocol-{}.jsonl'.format(rank)))
        current, active, setup = trace_state(events, rank, ranks)
        state = ('inside_exchange' if active else 'between_exchanges' if setup and current else
                 'setup_before_first_exchange' if setup else 'setup_returned')
        summary.update(records_valid=True, state=state, rejection='none')
        return summary, current if active else None
    except FileNotFoundError:
        summary['rejection'] = 'missing_file'
    except ValueError as error:
        # Only literal parser errors may cross the public boundary.
        summary['rejection'] = TRACE_REJECTIONS.get(str(error), 'invalid_record')
    except OSError:
        pass
    return summary, None


def report(directory, ranks):
    result = dict(assessment='incomplete_trace', failed_check='rank_scope', metadata_verified=False,
                  message_payloads_inspected=False, protocol_cause_verified=False,
                  scope='instrumented_host_message_calls_not_delivery_or_deadlock_proof')
    try:
        if type(ranks) is not int or not 1 <= ranks <= 4:
            raise ValueError('unsupported ranks')
        traces, summaries = [], {}
        for rank in range(ranks):
            summary, trace = rank_trace(directory, rank, ranks)
            summaries[str(rank)] = summary
            traces.append(trace)
        result['rank_trace_checks'] = summaries
        result['failed_check'] = 'trace_records'
        if not all(s['records_valid'] for s in summaries.values()):
            return result
        result['metadata_verified'] = True
        if any(t is None for t in traces):
            result.update(assessment='not_all_ranks_inside_exchange', failed_check='active_exchange')
            return result
        result.update(metadata_verified=True,
                      same_exchange_ordinal=len({t['entry']['exchange'] for t in traces}) == 1,
                      same_exchange_round=len({t['entry']['round'] for t in traces}) == 1,
                      all_communications_use_world=all(t['entry']['comm_world'] and all(
                          m['same_comm'] for m in t['messages'].values()) for t in traces),
                      mpi_call_error_seen=any(t['errors'] for t in traces))
        if not all(result[k] for k in ('same_exchange_ordinal', 'same_exchange_round', 'all_communications_use_world')):
            result.update(assessment='exchange_context_differs_or_unverified', failed_check='exchange_context')
            return result
        result['failed_check'] = 'wait_phase'
        round_number = traces[0]['entry']['round']
        sent, returned, received = Counter(), Counter(), Counter()
        for rank, trace in enumerate(traces):
            for ident, message in trace['messages'].items():
                key = (rank, message['peer'], message['tag'])
                if message['event'] in ('isend', 'send'):
                    sent[key] += 1
                    if trace['returned'].get(ident) is True:
                        returned[key] += 1
                elif message['event'] == 'recv' and trace['returned'].get(ident) is True:
                    received[(message['peer'], rank, message['tag'])] += 1
        result['rank_waits'] = {}
        for rank, trace in enumerate(traces):
            wait = trace['wait']
            blocking = any(i not in trace['returned'] for i in trace['messages'])
            phase, issued, calls_returned, contacts_received = 'unknown', None, None, None
            if wait and not blocking:
                requests = [trace['messages'][i] for i in wait['requests']]
                if requests and all(m['event'] == 'irecv' for m in requests):
                    tags = {m['tag'] for m in requests}
                    if tags == {1002 * round_number} and wait['event'] == 'testall':
                        phase = 'contact_replies'
                    elif tags == {1004 * round_number} and all(m['zero'] for m in requests):
                        phase = 'parent_release' if wait['event'] == 'test' else 'child_completion'
                    expected = Counter((m['peer'], rank, m['tag']) for m in requests)
                    issued = all(sent[k] >= n for k, n in expected.items())
                    calls_returned = all(returned[k] >= n for k, n in expected.items())
                    if phase == 'contact_replies':
                        contacts_received = all(received[(rank, peer, 1000 * round_number)] >= n
                                                for (peer, _, _), n in expected.items())
            result['rank_waits'][str(rank)] = dict(phase=phase,
                message_call_in_progress=blocking, matching_send_calls_seen=issued,
                matching_send_calls_returned=calls_returned,
                matching_contacts_received=contacts_received)
        result['assessment'] = 'protocol_waits_located' if all(
            r['phase'] != 'unknown' for r in result['rank_waits'].values()) else 'protocol_wait_unresolved'
        if result['assessment'] == 'protocol_waits_located':
            result['failed_check'] = 'none'
    except Exception:
        # Local paths, message partners, counters and raw exception text stay private.
        pass
    return result
