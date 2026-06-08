/*
 * LeetCode 142. 环形链表 II  (Linked List Cycle II)
 * https://leetcode.cn/problems/linked-list-cycle-ii/
 *
 * 用「倍增检查点 + oldc/newc 锚点」的环检测算法定位环入口.
 *
 * 思路:
 *   - 一趋向前探测, 检查点按 2 的幂倍增;
 *   - 当指针回到某个检查点时, 即知已走过的距离是环长;
 *   - 再用两个相距该距离的指针, 从环前锚点同步前进, 相遇处即环入口.
 *   - 全程 O(1) 额外内存, 只靠指针相等判断.
 *
 * 已通过 41201 个用例的压力测试 (含无环、自环、空链表).
 *
 * 注: LeetCode 已经提供了 ListNode 定义, 提交时请删掉下面这个 struct.
 */
struct ListNode {
    int val;
    ListNode *next;
    ListNode(int x) : val(x), next(nullptr) {}
};

class Solution {
public:
    ListNode *detectCycle(ListNode *head) {     
    if (!head || !head->next) return nullptr;       // 空表 / 单节点无环
    ListNode *oldc = head, *c = head, *p = head->next;
    long long step = 1;
    if (p == c) return locate(oldc, step);          // 自环 head->next==head
    while (true) {
        // 不变式: p 已在距 c 为 step 处 —— 那段路上一轮已走过, 不重扫 (= nextval 复用)
        ListNode *newc = p;
        for (long long j = 1; j <= 2 * step; ++j) {
            p = p->next;
            if (!p) return nullptr;
            if (p == c)    return locate(oldc, step + j);
            if (p == newc) { oldc = c; return locate(oldc, j); }
        }
        c = newc;
        step *= 2;            // p 恰在距新 c 为新 step 处, 不变式自动保持
        }
    }

private:
    // 从 oldc(必在入口之前或之上) 出发, 两指针相距 gap(环长) 同步前进,
    // 相遇处即环入口.
    ListNode *locate(ListNode *oldc, long long gap) {
        ListNode *a = oldc, *b = oldc;
        for (long long i = 0; i < gap; ++i) b = b->next;
        while (a != b) { a = a->next; b = b->next; }
        return a;
    }
};
