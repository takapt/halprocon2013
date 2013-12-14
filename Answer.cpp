//----------------------------------------------------------
/// @file
/// @brief    HPCAnswer.hpp の実装 (解答記述用ファイル)
/// @author   ハル研究所プログラミングコンテスト実行委員会
///
/// @copyright  Copyright (c) 2013 HAL Laboratory, Inc.
/// @attention  このファイルの利用は、同梱のREADMEにある
///             利用条件に従ってください

//----------------------------------------------------------

// Answer.cpp 専用のインクルードファイルです。
// 別のファイルをインクルードした場合、評価時には削除されます。
#include "HPCAnswerInclude.hpp"

#ifdef LOCAL
#include <bits/stdc++.h>
#include "local.h"
#else
#define assert(a)
#endif

// #define RUNTIME_DEBUG


#define rep(i, n) for (int i = 0; i < (int)(n); ++i)

typedef unsigned int uint;
typedef unsigned short ushort;
typedef unsigned char uchar;

namespace solver
{
using namespace hpc;

#ifdef LOCAL
using namespace std;
#endif



//// 
template <typename T>
void swap(T& a, T& b)
{
    T c = a;
    a = b;
    b = c;
}

template <typename T>
T min(const T& a, const T& b)
{
    return a < b ? a : b;
}

template <typename T>
T max(const T& a, const T& b)
{
    return a > b ? a : b;
}

template <typename T>
T abs(const T& n) { return n >= 0 ? n : -n; }

template <typename T>
void reverse(T* begin, T* end)
{
    --end;
    while (begin < end)
    {
        swap(*begin, *end);
        ++begin;
        --end;
    }
}

#ifdef RUNTIME_DEBUG
uint pss(uint first, uint second)
{
    assert(first == (ushort)(first));
    assert(second == (ushort)(second));
    return (first << 16) | second;
}
ushort pcc(uchar first, uchar second)
{
    assert(first == (uchar)(first));
    assert(second == (uchar)(second));
    return (first << 8) | second;
}

#else
uint pss(uint first, uint second) { return (first << 16) | second; }
ushort pcc(uchar first, uchar second) { return (first << 8) | second; }
#endif

ushort pss_first(uint pair) { return pair >> 16; }
ushort pss_second(uint pair) { return pair & 0xffff; }
uchar pcc_first(ushort pair) { return pair >> 8; }
uchar pcc_second(ushort pair) { return pair & 0xff; }



template <typename T, int SIZE>
class PriorityQueue
{
public:
    PriorityQueue()
        : n(0)
    {
    }

    void push(const T& a)
    {
        data[n++] = a;
        up(n - 1);
        assert(n <= SIZE);
    }

    void pop()
    {
        data[0] = data[--n];
        down(0);
    }

    T top() const
    {
        return data[0];
    }

    int size() const
    {
        return n;
    }

    bool empty() const
    {
        return size() == 0;
    }

    void clear()
    {
        n = 0;
    }

private:
    T data[SIZE];
    int n;

    void up(int k)
    {
        while (k > 0)
        {
            int par = (k - 1) / 2;
            if (!(data[k] < data[par]))
                break;

            swap(data[k], data[par]);
            k = par;
        }
    }

    void down(int k)
    {
        for (;;)
        {
            int next = k;
            int a = 2 * k + 1, b = 2 * k + 2;
            if (a < n && data[a] < data[next])
                next = a;
            if (b < n && data[b] < data[next])
                next = b;

            if (next == k)
                break;

            swap(data[k], data[next]);
            k = next;
        }
    }
};

class Random
{
private:
    unsigned int  x, y, z, w;
public:
    Random(unsigned int _x
             , unsigned int _y
             , unsigned int _z
             , unsigned int _w)
        : x(_x), y(_y), z(_z), w(_w) { }
    Random() 
        : x(123456789), y(362436069), z(521288629), w(88675123) { }
    Random(unsigned int seed)
        : x(123456789), y(362436069), z(521288629), w(seed) { }

    unsigned int next()
    {
        unsigned int t = x ^ (x << 11);
        x = y;
        y = z;
        z = w;
        return w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));
    }

    int next_int() { return next(); }

    // [0, upper)
    int next_int(int upper) { return next() % upper; }

    // [low, high]
    int next_int(int low, int high) { return next_int(high - low + 1) + low; }

    double next_double(double upper) { return upper * next() / 0xffffffffu; }
    double next_double(double low, double high) { return next_double(high - low) + low; }
};


/// 静止している2つの円が交差しているかどうかを返します。
/// 接している場合も交差とみなします。
///
/// @param[in] aC0    1つ目の円
/// @param[in] aC1    2つ目の円
///
/// @return 2つの円が衝突していたら @c true 。そうでなければ @c false 。
///         接しているときも衝突とみなします。
bool IsHit(const Circle& aC0, const Circle& aC1)
{
    const float squareDist = aC0.pos().squareDist(aC1.pos());
    return squareDist <= (aC0.radius() + aC1.radius()) * (aC0.radius() + aC1.radius());
}

/// 静止している円と、動いている円との衝突を判定します。
/// 接している場合も交差とみなします。
///
/// @param[in] aC0    静止している円。
/// @param[in] aC1    動いている円。移動開始時の状態を設定します。
/// @param[in] aP1    動いている円の移動後の位置。
///
/// @return 2つの円が衝突していたら @c true 。そうでなければ @c false 。
///         接しているときも衝突とみなします。
bool IsHit(const Circle& aC0, const Circle& aC1, const Vec2& aP1)
{
    // 動いていなければ静止している円の判定
    if(aC1.pos() == aP1) {
        return IsHit(aC0, aC1);
    }

    // 実質的に半径 aC0.radius + aC1.radius の円と線分 [aC1.pos, aP1.pos] の衝突判定と
    // 同等になる。
    const Circle c(aC0.pos(), aC0.radius() + aC1.radius());
    const Vec2 seg = aP1 - aC1.pos();
    const Vec2 c1ToC0 = c.pos() - aC1.pos();
    const float dist = Math::Abs(seg.cross(c1ToC0)) / seg.length();
    // 距離が c.radius より遠ければ衝突しない
    if(dist > c.radius()) {
        return false;
    }
    // 線分の延長線上で交差していないか調べる。
    const Vec2 p1ToC0 = c.pos() - aP1;
    // それぞれの点が円の反対方向にあれば衝突
    if(c1ToC0.dot(seg) * p1ToC0.dot(seg) <= 0.0f) {
        return true;
    }
    // 半径が大きければ衝突
    if(c.radius() >= c1ToC0.length() || c.radius() >= p1ToC0.length()) {
        return true;
    }
    return false;
}


class Action
{
public:
    Action(ActionType aType, float aValue)
        : mType(aType), mValue(aValue) {}
    Action()
        : mType(ActionType_TERM), mValue(-1) {}
    Action(const hpc::Action& act)
        : mType(act.type()), mValue(act.value()) {}

    ActionType type() const { return mType; }
    float value() const { return mValue; }

private:
    ActionType mType;
    float mValue;
};

inline bool valid_angle(float angle)
{
    return -Math::PI <= angle && angle <= Math::PI;
}
static int na = 0;
inline float normalize_angle(float angle)
{
    while (angle > Math::PI)
        angle -= 2 * Math::PI, ++na;
    while (angle < -Math::PI)
        angle += 2 * Math::PI, ++na;
    return angle;
}
inline float get_angle(const Vec2& pos)
{
    float ang = Math::ATan2(pos.y, pos.x);
    assert(valid_angle(ang));
    return ang;
}
inline float get_angle(float angle, float dest_angle)
{
    return normalize_angle(dest_angle - angle);
}
inline float get_angle(const Vec2& pos, float angle, const Vec2& dest)
{
    float dest_angle = get_angle(dest - pos);
    assert(valid_angle(dest_angle));
    return normalize_angle(get_angle(angle, dest_angle));
}

inline int need_move(float dist)
{
    assert(dist >= 0);
    const float eps = 1e-7;
    if (dist < eps)
        return 0;
    else
    {
        // return max((int)(Math::Ceil(dist / Parameter::PLAYER_SPEED - eps)), 1);
        return (int)(dist) + 1;
    }
}
inline int need_rot(float angle)
{
    angle = abs(angle);
    const float eps = 1e-7;
    if (angle < eps)
        return 0;
    else
    {
        // return max((int)(Math::Ceil(angle / Parameter::PLAYER_ANGLE_RATE - eps)), 1);
        return (int)(angle / Parameter::PLAYER_ANGLE_RATE) + 1;
    }
}

// class Player
// {
// public:
//     Player()
//         : _region(Circle(Vec2(1e4, 1e4), 1e4)), _angle(1e4)
//     {
//     }
//     explicit
//     Player(const hpc::Player& player)
//         : _region(player.region()), _angle(normalize_angle(player.arg()))
// #ifdef RUNTIME_DEBUG
//           ,_player(player)
// #endif
//     {
//         _region.setRadius(Parameter::PLAYER_RADIUS - 1e-5);
//     }
// 
//     void act(const Action& action)
//     {
// #ifdef RUNTIME_DEBUG
//         _player.act(hpc::Action(action.type(), action.value()));
// #endif
//         switch (action.type())
//         {
//             case ActionType_Move:
//                 move(action.value());
//                 break;
//             case ActionType_Rotate:
//                 rotate(action.value());
//                 break;
//             default:
//                 assert(false);
//         }
// 
// #ifdef RUNTIME_DEBUG
//         // if (!Math::IsEqual(arg(), normalize_angle(_player.arg())))
//         if (abs(normalize_angle(arg() - normalize_angle(_player.arg()))) > 1e-4)
//         {
//             printf("not eq: %.10f, %.10f\n", arg(), normalize_angle(_player.arg()));
//             exit(1);
//         }
//         // if (!pos().equals(_player.pos()))
//         if (abs(pos().x - _player.pos().x) > 1e-2 ||
//             abs(pos().y - _player.pos().y) > 1e-2)
//         {
//             printf("not eq: (%f, %f), (%f, %f)\n", pos().x, pos().y,
//                     _player.pos().x, _player.pos().y);
//             exit(1);
//         }
// #endif
//     }
// 
//     void move(float dist)
//     {
//         const float dx = dist * Math::Cos(arg());
//         const float dy = dist * Math::Sin(arg());
//         _region.move(Vec2(dx, dy));
//     }
// 
//     void rotate(float rot)
//     {
//         assert(valid_angle(rot)); // for fast
//         _angle += rot;
//         // _angle = normalize_angle(_angle);
//         if (_angle > Math::PI)
//             _angle -= 2 * Math::PI;
//         else if (_angle < -Math::PI)
//             _angle += 2 * Math::PI;
//         assert(valid_angle(_angle));
//     }
// 
//     const Circle& region() const { return _region; }
// 
//     Vec2 vec() const
//     {
//         return Vec2(Math::Cos(arg()), Math::Sin(arg()));
//     }
// 
//     Vec2 pos() const { return _region.pos(); }
// 
//     float arg() const
//     {
//         // if (!valid_angle(_angle))
//         // {
//         //     debug(_angle);
//         // }
//         // assert(valid_angle(_angle));
//         return _angle;
//     }
// 
// private:
//     Circle _region;
//     float _angle;
// 
// #ifdef RUNTIME_DEBUG
//     hpc::Player _player;
// #endif
// };

class Answer
{
public:
    Answer()
        : _size(0)
    {
    }

    void reset() { _size = 0; }

    int size() const { return _size; }
    Action action(int i) const
    {
        assert(0 <= i && i < size());
        return _actions[i];
    }

    void add(const Action& act)
    {
        assert(_size < Parameter::GAME_TURN_PER_STAGE);
        _actions[_size++] = act;
    }

#ifdef RUNTIME_DEBUG
    void add(const Action& act, const solver::Player& p)
    {
        _player_log[_size] = p;
        add(act);
    }
#endif

    int num_move() const
    {
        int res = 0;
        rep(i, size())
            if (action(i).type() == ActionType_Move)
                ++res;
        return res;
    }
    int num_rot() const
    {
        int res = 0;
        rep(i, size())
            if (action(i).type() == ActionType_Rotate)
                ++res;
        return res;
    }

private:
    int _size;
    Action _actions[Parameter::GAME_TURN_PER_STAGE];

#ifdef RUNTIME_DEBUG
public:
    solver::Player player(int i) const { return _player_log[i]; }
private:
    solver::Player _player_log[Parameter::GAME_TURN_PER_STAGE];
#endif
};



inline bool intersect(const Vec2& a, const Vec2& b, const Vec2& c, const Vec2& d)
{
    Vec2 ab = b - a;
    if (ab.cross(c - a) * ab.cross(d - a) > 0)
        return false;
    Vec2 cd = d - c;
    if (cd.cross(a - c) * cd.cross(b - c) > 0)
        return false;
    return true;
}
inline bool can_go_straight(const Vec2& a, const Vec2& b, const HoleCollection& holes)
{
    rep(i, holes.count())
    {
        if (intersect(a, b, holes[i].pointA(), holes[i].pointB()))
            return false;
        if (intersect(a, b, holes[i].pointB(), holes[i].pointC()))
            return false;
        if (intersect(a, b, holes[i].pointC(), holes[i].pointD()))
            return false;
        if (intersect(a, b, holes[i].pointD(), holes[i].pointA()))
            return false;
    }
    return true;
}


Rectangle expand_margin(const Rectangle& rect, float margin)
{
    return Rectangle(rect.center(), rect.width() + 2 * margin, rect.height() + 2 * margin);
}

int stage_no = -1;
int calc_cost_call = 0;
class OrderSolver
{
public:
    OrderSolver(const Stage& _stage)
        : stage(_stage), items(stage.items()), n(items.count()), player(_stage.player())
    {
        rep(i, n)
            order[i] = i;
    }

    void get_order(int* res) const
    {
        rep(i, n)
            res[i] = order[i];
    }

    void solve()
    {
        greedy_near();
        improve();
    }

private:
    void greedy_near()
    {
        bool used[128] = {};

        Vec2 pos = stage.player().pos();
        rep(oi, n)
        {
            int k = -1;
            float min_d = 1e7;
            rep(i, n)
            {
                if (!used[i])
                {
                    float d = pos.dist(items[i].pos());
                    if (d < min_d)
                    {
                        min_d = d;
                        k = i;
                    }
                }
            }
            assert(k != -1);

            order[oi] = k;
            pos = items[k].pos();
            used[k] = true;
        }
    }

    void improve()
    {
        if (n == 1)
            return;

        static Random rand;
        int cur = calc_cost();
        // for (int loop = 0; loop < n * n * 10; ++loop)
            // int a = rand.next_int(n);
            // int b = (a + rand.next_int(n - 1)) % n;
            // swap(order[i], order[j]);

            // insert
            // bool updated = false;
            // rep(b, n) rep(a, b)
            // {
            //     int ori[128];
            //     rep(i, n)
            //         ori[i] = order[i];
            //     if (a > b)
            //         swap(a, b);
            //     rep(i, n)
            //         order[i] = ori[i];
            //     for (int i = a; i < min(n - 1, b); ++i)
            //         swap(order[i], order[i + 1]);

            //     int next = calc_cost();
            //     if (next < cur)
            //     {
            //         cur = next;
            //         updated = true;
            //     }
            //     else
            //     {
            //         // swap(order[i], order[j]);
            //         rep(i, n)
            //             order[i] = ori[i];
            //     }
            // }
            // if (!updated)
            //     break;

        for (bool updated = true; updated; )
        {
            updated = false;
            rep(b, n - 1) rep(a, b)
            {
                int ori[128];
                rep(i, n)
                    ori[i] = order[i];

                swap(order[a + 1], order[b]);
                reverse(order + a + 2, order + b);

                int next = calc_cost();
                if (next < cur)
                {
                    cur = next;
                    updated = true;
                }
                else
                {
                    // swap(order[i], order[j]);
                    rep(i, n)
                        order[i] = ori[i];
                }
            }
        }
    }

    // k周りにかかるコスト
    // int cost_around(int k)
    // {
    //     // move
    //     int move_cost = 0;
    //     if (k > 0)
    //         move_cost += need_move(items[order[k - 1]].pos().dist(items[order[k]].pos()));
    //     if (k < n - 1)
    //         move_cost += need_move(items[order[k]].pos().dist(items[order[k + 1]].pos()));

    //     // rotation
    //     int num_p = 0;
    //     Vec2 points[5];
    //     float angle;
    //     if (k == 0)
    //     {
    //         points[num_p++] = player.pos();
    //         angle = player.arg();
    //     }
    //     else if (k == 1)
    //     {
    //         points[num_p++] = items[order[0]].pos();
    //         angle = get_angle(items[order[0]].pos() - player.pos());
    //     }
    //     else
    //     {
    //         points[num_p++] = items[order[k - 1]].pos();
    //         angle = get_angle(items[order[k - 1]].pos() - items[order[k - 2]].pos());
    //     }
    //     points[num_p++] = items[order[k]].pos();
    //     if (k < n - 1)
    //         points[num_p++] = items[order[k + 1]].pos();

    //     int rot_cost = 0;
    //     rep(i, num_p - 1)
    //     {
    //         float next_angle = get_angle(points[i + 1] - points[i]);
    //         rot_cost += need_rot(get_angle(angle, next_angle));
    //         angle = next_angle;
    //     }

    //     int cost = move_cost + rot_cost;
    //     return cost;
    // }

    int calc_cost() const
    {
        ++calc_cost_call;

        const int num_p = n + 1;
        Vec2 points[128];
        points[0] = player.pos();
        rep(i, n)
            points[i + 1] = items[order[i]].pos();

        int move_cost = 0;
        rep(i, num_p - 1)
            move_cost += need_move(points[i].dist(points[i + 1]));

        int rot_cost = 0;
        rot_cost += need_rot(get_angle(player.pos(), player.arg(), items[order[0]].pos()));
        float angle = get_angle(points[1] - points[0]);
        for (int i = 1; i < num_p - 1; ++i)
        {
            float next_angle = get_angle(points[i + 1] - points[i]);
            rot_cost += need_rot(get_angle(angle, next_angle));
            angle = next_angle;
        }

        int cost = move_cost + rot_cost;
        return cost;
    }

private:
    const Stage& stage;
    const ItemCollection& items;
    const int n;
    const solver::Player player;

    int order[128];
};
void solve_order(const Stage& stage, int* order)
{
    OrderSolver s(stage);
    s.solve();
    s.get_order(order);
}

void solve_actions(const Stage& stage, int* item_order, Answer& answer)
{
    answer.reset();

    ItemCollection items = stage.items();
    const HoleCollection& holes = stage.holes();
    const int n = items.count();


    const float margin = 0.01;
    int num_corners = 0;
    Vec2 points[64];
    rep(i, holes.count())
    {
        Rectangle rect = expand_margin(holes[i], margin);
        points[num_corners++] = rect.pointA();
        points[num_corners++] = rect.pointB();
        points[num_corners++] = rect.pointC();
        points[num_corners++] = rect.pointD();
    }

    const float inf = 1e6;
    float dist[64][64];
    rep(i, num_corners) rep(j, num_corners)
        dist[i][j] = inf;

    rep(i, num_corners) rep(j, i)
    {
        if (can_go_straight(points[i], points[j], holes))
        {
            dist[i][j] = dist[j][i] = points[i].dist(points[j]);
        }
    }

    // hpc::Player player = stage.player();
    solver::Player player(stage.player());
    rep(oi, n)
    {
        // debug(oi);
        const Item& item = items[item_order[oi]];
        if (item.isRemoved())
            continue;

        int num_mark;
        Vec2 mark[64];
        if (can_go_straight(player.pos(), item.pos(), holes))
        {
            num_mark = 1;
            mark[0] = item.pos();
        }
        else
        {
            const int start = num_corners;
            const int goal = num_corners + 1;
            int num_points = num_corners;
            points[num_points++] = player.pos();
            points[num_points++] = item.pos();
            rep(i, num_points)
            {
                for (int j = num_corners; j < num_points; ++j)
                {
                    if (i != j && can_go_straight(points[i], points[j], holes))
                        dist[i][j] = dist[j][i] = points[i].dist(points[j]);
                    else
                        dist[i][j] = dist[j][i] = inf;
                }
            }

            // 順番を求める
            const int INT_INF = inf;
            int dp[64][64]; // (pos, next_dest_pos)
            int prev[64][64];
            rep(i, num_points) rep(j, num_points)
            {
                dp[i][j] = INT_INF;
                prev[i][j] = -INT_INF;
            }
            PriorityQueue<uint, 64 * 64> q;
            rep(i, num_points)
            {
                if (dist[start][i] < inf)
                {
                    int rot_cost = need_rot(get_angle(player.pos(), player.arg(), points[i]));
                    dp[start][i] = rot_cost;
                    q.push(pss(rot_cost, pcc(start, i)));
                }
            }

            while (!q.empty())
            {
                int cost = pss_first(q.top());
                int pos = pcc_first(pss_second(q.top()));
                int npos = pcc_second(pss_second(q.top()));
                q.pop();

                // printf("%d, %d: %d\n", pos, npos, cost);
                assert(0 <= pos && pos < num_points);
                assert(0 <= npos && npos < num_points);

                if (cost > dp[pos][npos])
                    continue;

                if (pos == goal)
                {
                    break;
                }

                cost += need_move(dist[pos][npos]);
                if (npos == goal)
                {
                    if (cost < dp[npos][0])
                    {
                        dp[npos][0] = cost;
                        prev[npos][0] = pos;
                        q.push(pss(cost, pcc(npos, 0)));
                    }
                }
                else
                {
                    float cur_angle = get_angle(points[npos] - points[pos]);
                    rep(i, num_points)
                    {
                        if (dist[npos][i] < inf)
                        {
                            int ncost = cost + need_rot(get_angle(cur_angle, get_angle(points[i] - points[npos])));
                            if (ncost < dp[npos][i])
                            {
                                dp[npos][i] = ncost;
                                prev[npos][i] = pos;
                                q.push(pss(ncost, pcc(npos, i)));
                            }
                        }
                    }
                }
            }
            assert(dp[goal][0] != INT_INF);

            int num = 0;
            int order[64];
            for (int pos = goal, npos = 0; pos != start; )
            {
                order[num++] = pos;
                assert(num <= num_points);

                int _pos = pos;
                pos = prev[pos][npos];
                npos = _pos;
            }
            assert(num <= num_points);
            reverse(order, order + num);

            num_mark = num;
            rep(i, num)
                mark[i] = points[order[i]];
            if (oi + 1 < n)
                mark[num_mark] = items[oi + 1].pos();
        }

        // action
        rep(mi, num_mark)
        {
            const Vec2& dest = mark[mi];

            while (!item.isRemoved() && !player.pos().equals(dest))
            {
                // printf("(%.3f, %.3f, %.3f), (%.3f, %.3f)\n",
                //         player.pos().x, player.pos().y, Math::RadToDeg(player.arg()), dest.x, dest.y);

                const float eps = 1e-4;

                const Vec2& cur = player.pos();
                const float& cur_angle = normalize_angle(player.arg());
                const float d = cur.dist(dest);
                const float angle_diff = get_angle(cur, cur_angle, dest);
                float dest_angle = angle_diff;

                bool force_move = false;
                const Vec2& cur_vec = player.vec();
                Item titem = item;
                if (mi == num_mark - 1 && oi < n - 1 &&
                    // IsHit(item.region(), player.region(), cur + cur_vec) &&
                    titem.tryToRemove(player.region(), 100 * cur_vec) &&
                    cur_vec.cross(dest - cur) * cur_vec.cross(items[oi + 1].pos() - cur) < 0 &&
                    stage.field().isIn(cur + d * cur_vec) &&
                    can_go_straight(cur, cur + d * cur_vec, holes)
                   )
                {
                    assert(dest.equals(item.pos()));
                    force_move = true;
                    // cout << "force" << endl;
                }

                if (!force_move && abs(dest_angle) > eps)
                {
                    float rot = Math::LimitAbs(dest_angle, Parameter::PLAYER_ANGLE_RATE);
                    hpc::Action action(ActionType_Rotate, rot);

#ifdef RUNTIME_DEBUG
                    answer.add(action, player);
#else
                    answer.add(action);
#endif

                    player.act(action);
                }
                else
                {
                    float move = min(d, Parameter::PLAYER_SPEED);
                    hpc::Action action(ActionType_Move, move);

                    rep(i, items.count())
                    {
                        items[i].tryToRemove(player.region(), player.vec() * move);
                        // items[i].tryToRemove(Circle(player.pos(), Parameter::PLAYER_RADIUS), player.vec() * move);
                    }

#ifdef RUNTIME_DEBUG
                    answer.add(action, player);
#else
                    answer.add(action);
#endif

                    player.act(action);

                }
            }
            // cout << "ouou" << endl;
        }
        // cout << "loopend" << endl;
    }
}


Answer answer;
int ans_i;
void solve(const Stage& stage)
{
    int order[128];
    solve_order(stage, order);
    solve_actions(stage, order, answer);

    ans_i = 0;
}

hpc::Action get_next_action()
{
    const Action& t = answer.action(ans_i++);
    return hpc::Action(t.type(), t.value());
}

#ifdef RUNTIME_DEBUG
Player get_player()
{
    return answer.player(ans_i);
}
#endif
}


/// プロコン問題環境を表します。
namespace hpc {

    //----------------------------------------------------------
    /// 各ステージ開始時に呼び出されます。
    ///
    /// この関数を実装することで、各ステージに対して初期設定を行うことができます。
    ///
    /// @param[in] aStage 現在のステージ。
    void Answer::Init(const Stage& aStage)
    {
        ++solver::stage_no;
        debug(solver::stage_no);
        solver::solve(aStage);
        // printf("%d, %d\n", solver::answer.num_move(), solver::answer.num_rot());
        // debug(solver::na);
        // debug(solver::calc_cost_call);
        static int move = 0, rot = 0;
        move += solver::answer.num_move();
        rot += solver::answer.num_rot();
        // printf("%4d %4d (%.3f) %8d %8d\n",
        //         solver::answer.num_move(), solver::answer.num_rot(),
        //         float(solver::answer.num_rot()) / solver::answer.num_move(),
        //         move, rot);
    }

    //----------------------------------------------------------
    /// 各ターンでの動作を返します。
    ///
    /// @param[in] aStage 現在ステージの情報。
    ///
    /// @return これから行う動作を表す Action クラス。
    ///
    /// @attention 戻り値の Action の値が、それぞれの動作量の最大値・最小値を外れる場合、
    ///            自動的に最大値または最小値に制限されます。
    Action Answer::GetNextAction(const Stage& aStage)
    {
        // Player p = aStage.player();
        // Player sp = solver::get_player();
        // float dx = sp.pos().x - p.pos().x;
        // float dy = sp.pos().y - p.pos().y;
        // float da = solver::normalize_angle(sp.arg() - p.arg());
        // if (abs(dx) > 1e-4 || abs(dy) > 1e-4 || abs(da) > 1e-4)
        // {
        //     printf("(%3.3f, %3.3f, %3.3f)\t(%3.3f, %3.3f, %3.3f)\n",
        //             p.pos().x, p.pos().y, p.arg(),
        //             sp.pos().x, sp.pos().y, sp.arg());
        //     printf("error: %3.10f, %3.10f, %3.10f\n", dx, dy, da);
        //     puts("");
        // }
        return solver::get_next_action();
    }
}

//----------------------------------------------------------
// EOF
