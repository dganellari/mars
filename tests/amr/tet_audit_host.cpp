#include "backend/distributed/unstructured/amr/mars_amr_tet_audit_math.hpp"
#include "core/mars_tet_red.hpp"
#include <algorithm>
#include <array>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace mars::amr;

void require(bool condition, const char* message)
{
    if (!condition) throw std::runtime_error(message);
}

int main()
{
    double points[10][3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    int edge             = 4;
    for (int a = 0; a < 4; ++a)
        for (int b = a + 1; b < 4; ++b, ++edge)
            for (int d = 0; d < 3; ++d)
                points[edge][d] = (points[a][d] + points[b][d]) / 2;
    uint64_t nodes[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    std::vector<TetAuditRecord<uint64_t>> children;
    TetAuditTotals total;
    for (int child = 0; child < 8; ++child)
    {
        uint64_t keys[4];
        red_child_connectivity(child, nodes, keys);
        double p[4][3];
        for (int i = 0; i < 4; ++i)
            for (int d = 0; d < 3; ++d)
                p[i][d] = points[keys[i]][d];
        auto record = audit_record(keys, p);
        require(!record.invalid && record.volume > 0, "red child validity/orientation");
        children.push_back(record);
        total = TetAuditFold{}(total, TetAuditValue<uint64_t>{}(record));
    }
    require(total.count == 8 && std::abs(total.volume - 1.0 / 6) < 1e-15, "red volume conservation");
    TetAuditTotals reversed;
    for (auto i = children.rbegin(); i != children.rend(); ++i)
        reversed = TetAuditFold{}(reversed, TetAuditValue<uint64_t>{}(*i));
    require(audit_same_cells(total, reversed), "redistribution/order independence");
    TetAuditTotals missing, replaced;
    for (size_t i = 0; i < children.size(); ++i)
    {
        if (i) missing = TetAuditFold{}(missing, TetAuditValue<uint64_t>{}(children[i]));
        replaced = TetAuditFold{}(replaced, TetAuditValue<uint64_t>{}(children[i ? i : 1]));
    }
    require(!audit_same_cells(total, missing), "missing cell must fail");
    require(replaced.count == total.count && std::abs(replaced.volume - total.volume) < 1e-15,
            "duplicate/loss control preserves counts and volume");
    require(!audit_same_cells(total, replaced), "duplicate/loss must fail identity check");

    uint64_t keys[4] = {9, 2, 7, 3};
    double p[4][3]   = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    auto positive    = audit_record(keys, p);
    require(audit_red_children_valid(positive, children.data(), children.size(), 0, true), "per-parent children");
    auto damaged = children;
    damaged[0] = damaged[1];
    require(!audit_red_children_valid(positive, damaged.data(), damaged.size(), 0, true),
            "duplicated child with unchanged volume must fail");
    damaged = children;
    damaged[0].volume *= 1.1;
    damaged[1].volume -= damaged[0].volume - children[0].volume;
    require(!audit_red_children_valid(positive, damaged.data(), damaged.size(), 0, true),
            "compensating child volume errors must fail");
    require(!audit_red_children_valid(positive, children.data(), children.size() - 1, 0, true),
            "incomplete child range must fail");
    std::swap(keys[0], keys[1]);
    for (int d = 0; d < 3; ++d)
        std::swap(p[0][d], p[1][d]);
    auto negative = audit_record(keys, p);
    require(audit_same_cells(TetAuditValue<uint64_t>{}(positive), TetAuditValue<uint64_t>{}(negative)),
            "canonical identity ignores corner ordering");
    require(TetAuditValue<uint64_t>{}(negative).inverted == 1 && negative.volume == -positive.volume,
            "orientation remains separately observable");
    keys[1] = keys[0];
    require(audit_record(keys, p).invalid, "duplicate key must fail");
    keys[1] = 4;
    p[0][0] = std::numeric_limits<double>::quiet_NaN();
    require(audit_record(keys, p).invalid, "NaN must fail");
    p[0][0] = 0;
    for (auto& point : p)
        point[2] = 0;
    require(audit_record(keys, p).invalid, "collapsed tet must fail");
    std::cout << "PASS: red children, volume, ownership order, omissions, compensating duplicates, orientation and "
                 "invalid records\n";
}
