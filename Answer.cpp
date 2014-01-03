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

namespace std
{
ostream& operator<<(ostream& os, const hpc::Vec2& p)
{
    char buf[128];
    sprintf(buf, "(%7.3f, %7.3f)", p.x, p.y);
    os << buf;
    return os;
}
}
#else
#define assert(a)
#define abort()
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

const char* to_s(const Vec2& p)
{
    static int i = 0;
    ++i %= 1000;
    static char buf[1024][128];
    sprintf(buf[i], "(%3.3f, %3.3f)", p.x, p.y);
    return buf[i];
}
#endif



//// 
template <typename T>
inline
void swap(T& a, T& b)
{
    T c = a;
    a = b;
    b = c;
}

template <typename T>
inline
T min(const T& a, const T& b)
{
    return a < b ? a : b;
}

template <typename T>
inline
T max(const T& a, const T& b)
{
    return a > b ? a : b;
}

template <typename T>
inline
T abs(const T& n) { return n >= 0 ? n : -n; }

template <typename T>
inline
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

////////////////////////////////// for AI
int stage_no = -1;

// 点pから直線abへの垂線の足
Vec2 projection(const Vec2& p, const Vec2& a, const Vec2& b)
{
    const Vec2 b_to_a = a - b;
    float t = (p - a).dot(b_to_a) / b_to_a.squareLength();
    return a + t * b_to_a;
}
// 円cと直線abとのa側の交点
Vec2 intersect_point(const Circle& c, const Vec2& a, const Vec2& b)
{
    Vec2 p = projection(c.pos(), a, b);
    Vec2 c_to_p = p - c.pos();
    float d = c_to_p.length();
    assert(d < c.radius() + 1e-3); // 円と直線は交わってる
    if (d > c.radius()) // assertに引っかからない範囲なら接線との交点を返しとく
        return p;

    float t = Math::Sqrt(c.radius() * c.radius() - d * d);
    const Vec2 u = b - a;
    const Vec2 v = t / u.length() * u;
    return p - v;
}
void intersect_points(const Circle& c, const Vec2& a, const Vec2& b, Vec2& ip1, Vec2& ip2)
{
    Vec2 p = projection(c.pos(), a, b);
    Vec2 c_to_p = p - c.pos();
    float d = c_to_p.length();
    assert(d < c.radius() + 1e-3); // 円と直線は交わってる
    if (d > c.radius()) // assertに引っかからない範囲なら接線との交点を返しとく
    {
        ip1 = ip2 = p;
        return;
    }

    float t = Math::Sqrt(c.radius() * c.radius() - d * d);
    const Vec2 u = b - a;
    const Vec2 v = t / u.length() * u;
    ip1 = p - v;
    ip2 = p + v;
}

bool intersect(const Circle& aC0, const Circle& aC1)
{
    const float squareDist = aC0.pos().squareDist(aC1.pos());
    return squareDist <= (aC0.radius() + aC1.radius()) * (aC0.radius() + aC1.radius());
}
bool intersect(const Circle& c, const Vec2& p)
{
    return p.squareDist(c.pos()) <= c.radius() * c.radius();
}
bool intersect(const Circle& c, const Vec2& from, const Vec2& to)
{
    // 動いていなければ静止している円の判定
    if (from == to)
        return intersect(c, from);

    // 実質的に半径 aC0.radius + aC1.radius の円と線分 [aC1.pos, aP1.pos] の衝突判定と
    // 同等になる。
    // const Circle c(aC0.pos(), aC0.radius() + aC1.radius());
    // const Vec2 seg = aP1 - aC1.pos();
    // const Vec2 c1ToC0 = c.pos() - aC1.pos();
    const Vec2 seg = to - from;
    const Vec2 from_to_c = c.pos() - from;
    // const float dist = Math::Abs(seg.cross(c1ToC0)) / seg.length();
    const float dist = Math::Abs(seg.cross(from_to_c)) / seg.length();
    // 距離が c.radius より遠ければ衝突しない
    if(dist > c.radius()) {
        return false;
    }
    // 線分の延長線上で交差していないか調べる。
    // const Vec2 p1ToC0 = c.pos() - aP1;
    const Vec2 to_to_c = c.pos() - to;
    // それぞれの点が円の反対方向にあれば衝突
    // if(c1ToC0.dot(seg) * p1ToC0.dot(seg) <= 0.0f) {
    if (from_to_c.dot(seg) * to_to_c.dot(seg) <= 0)
        return true;
    // 半径が大きければ衝突
    // if(c.radius() >= c1ToC0.length() || c.radius() >= p1ToC0.length()) {
    if (c.radius() >= to_to_c.length() || c.radius() >= to_to_c.length())
        return true;
    return false;
}
bool intersect(const Circle& aC0, const Circle& aC1, const Vec2& aP1)
{
    // 動いていなければ静止している円の判定
    if(aC1.pos() == aP1) {
        return intersect(aC0, aC1);
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

bool contain(const Rectangle& rect, const Circle& c)
{
    const Vec2 corner[] = { rect.pointA(), rect.pointB(), rect.pointC(), rect.pointD() };
    rep(i, 4)
    {
        if (intersect(c, corner[i], corner[(i + 1) % 4]))
            return false;
    }
    return true;
}

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


    // こっちのがわずかに遅い？
    // if (angle > Math::PI)
    // {
    //     int k = (int)((angle - Math::PI) / (2 * Math::PI)) + 1;
    //     angle -= k * 2 * Math::PI;
    // }
    // else if (angle < -Math::PI)
    // {
    //     int k = (int)(-(angle + Math::PI) / (2 * Math::PI)) + 1;
    //     angle += k * 2 * Math::PI;
    // }

    assert(valid_angle(angle));
    return angle;
}
inline float get_angle(const Vec2& pos)
{
    float ang = Math::ATan2(pos.y, pos.x);
    assert(valid_angle(ang));
    return ang;
}
inline float get_angle(const Vec2& pos_from, const Vec2& pos_to)
{
    return get_angle(pos_to - pos_from);
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
    angle = normalize_angle(angle);
    assert(valid_angle(angle));

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

// pos -> pos[0] -> pos[1] -> ... -> pos[n - 1]のコスト
int calc_cost(Vec2 cur, float angle, const Vec2* pos, int n)
{
    int total = 0;
    rep(i, n)
    {
        float dest_angle = get_angle(cur, pos[i]);
        total += need_rot(dest_angle - angle);
        total += cur.dist(pos[i]);

        cur = pos[i];
        angle = dest_angle;
    }
    return total;
}
// cur -> a -> bと行動するときのコスト
int calc_cost(const Vec2& cur, float angle, const Vec2& a, const Vec2& b)
{
    Vec2 pos[] = { a, b };
    return calc_cost(cur, angle, pos, 2);
}

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



bool intersect(const Vec2& a, const Vec2& b, const Vec2& c, const Vec2& d)
{
    Vec2 ab = b - a;
    if (ab.cross(c - a) * ab.cross(d - a) > 0)
        return false;
    Vec2 cd = d - c;
    if (cd.cross(a - c) * cd.cross(b - c) > 0)
        return false;
    return true;
}
bool intersect(const Vec2& a, const Vec2& b, const Rectangle& rect)
{
    if (intersect(a, b, rect.pointA(), rect.pointB()))
        return true;
    if (intersect(a, b, rect.pointB(), rect.pointC()))
        return true;
    if (intersect(a, b, rect.pointC(), rect.pointD()))
        return true;
    if (intersect(a, b, rect.pointD(), rect.pointA()))
        return true;
    return false;
}
bool can_go_straight(const Vec2& a, const Vec2& b, const HoleCollection& holes)
{
    rep(i, holes.count())
        if (intersect(a, b, holes[i]))
            return false;
    return true;
}


Rectangle expand_margin(const Rectangle& rect, float margin)
{
    return Rectangle(rect.center(), rect.width() + 2 * margin, rect.height() + 2 * margin);
}

// complex product
Vec2 product(const Vec2& a, const Vec2& b)
{
    return Vec2(a.x * b.x - a.y * b.y, a.y * b.x + a.x * b.y);
}

// lower.cross(upper) > 0
void tangent_vecs(const Vec2& point, const Circle& circle, Vec2& lower, Vec2& upper)
{
    Vec2 vec = circle.pos() - point;
    const float d = vec.length();
    assert(d > circle.radius()); // pointがcircle外

    const float r = circle.radius();
    const float t = Math::Sqrt(d * d - r * r);
    Vec2 rot(t / d, r / d);
    vec *= t / d; // adjust length

    upper = product(vec, rot);

    rot.y *= -1;
    lower = product(vec, rot);

    // cerr << endl;
    // debug(to_s(point));
    // printf("%s %f\n", to_s(circle.pos()).c_str(), circle.radius());
    assert(lower.cross(upper) > 0);
}

// return: hole iと交わる時、result |= 1 << i
int intersect_holes_mask(const Vec2& from, const Vec2& to, const HoleCollection& holes)
{
    int mask = 0;
    rep(i, holes.count())
    {
        if (intersect(from, to, holes[i]))
            mask |= 1 << i;
    }
    return mask;
}

// return: 有効な範囲があるか
bool and_range(const Vec2& lower_a, const Vec2& upper_a, const Vec2& lower_b, const Vec2& upper_b, Vec2& lower, Vec2& upper)
{
    if (lower_a.cross(upper_b) < 0 || upper_a.cross(lower_b) > 0)
        return false;

    lower = lower_a.cross(lower_b) > 0 ? lower_b : lower_a;
    upper = upper_a.cross(upper_b) < 0 ? upper_b : upper_a;

// #define RAN(low, up) Math::RadToDeg(get_angle(low)), Math::RadToDeg(get_angle(up))
//     printf("[%f, %f] and [%f, %f] -> [%f, %f]\n", RAN(lower_a, upper_a), RAN(lower_b, upper_b), RAN(lower, upper));
    return true;
}

class ActionSolver
{
public:
    ActionSolver(const Stage& _stage)
        : stage(_stage), holes(_stage.holes())
    {
        const float hole_margin = 0.001;
        rep(i, holes.count())
        {
            Rectangle rect = expand_margin(holes[i], hole_margin);
            corners[i][0] = rect.pointA();
            corners[i][1] = rect.pointB();
            corners[i][2] = rect.pointC();
            corners[i][3] = rect.pointD();
        }

        const float inf = 1e6;
        rep(i, holes.count()) rep(j, holes.count()) rep(a, 4) rep(b, 4)
            corner_move_cost[i][a][j][b] = inf;
        rep(j, holes.count()) rep(i, j + 1) rep(a, 4) rep(b, 4)
        {
            if ((i != j || a != b) && can_go_straight(corners[i][a], corners[j][b], holes))
            {
                corner_move_cost[i][a][j][b] = corner_move_cost[j][b][i][a] = need_move(corners[i][a].dist(corners[j][b]));
                corner_angle[i][a][j][b] = get_angle(corners[j][b] - corners[i][a]);
                corner_angle[j][b][i][a] = get_angle(corners[i][a] - corners[j][b]);
            }
        }

        const ItemCollection& _items = stage.items();
        rep(i, _items.count())
        {
            static const float eps = 1e-1;
            items.setupAddItem(_items[i].pos(), _items[i].radius() + Parameter::PLAYER_RADIUS - eps);
            ori_items.setupAddItem(_items[i].pos(), _items[i].radius() + Parameter::PLAYER_RADIUS);
        }
    }

public:
    int solve(int* item_order, Vec2* turning_points)
    {
        ItemCollection debug_items = stage.items();

        Vec2 cur_pos = stage.player().pos();
        float cur_angle = stage.player().arg();

        Vec2 next_dest[128]; // next_dest[i] -> items[i]を取った後に向かうべき座標
        {
            float angle = cur_angle;
            rep(i, items.count() - 1)
            {
                const Vec2& cur = items[item_order[i]].pos();
                const Vec2& next = items[item_order[i + 1]].pos();

                int holes_mask = intersect_holes_mask(cur, next, holes);
                if (holes_mask)
                    next_dest[item_order[i]] = next_dest_pos(cur, angle, next, holes_mask);
                else
                    next_dest[item_order[i]] = next;


                angle = get_angle(next_dest[item_order[i]], next);

                assert(can_go_straight(cur, next_dest[item_order[i]], holes));
            }
        }

        // cerr << endl << endl;
        // cerr << "item_order: ";
        // print(item_order, items.count());
        // rep(i, items.count())
        // {
        //     int k = item_order[i];
        //     fprintf(stderr, "%2d: ", k);
        //     cerr << next_dest[item_order[i]] << endl;
        // }

        int num_turning_points = 0;
        turning_points[num_turning_points++] = cur_pos;
        // for (int item_i = 0; item_i < items.count(); )
        for (int order_i = 0; order_i < items.count(); )
        {
            const int item_i = item_order[order_i];

            // cerr << endl;
            // cerr << string(30, '-') << endl;
            // debug(item_i);

            // assert(!intersect(items[item_i].region(), cur_pos));
            if (intersect(items[item_i].region(), cur_pos))
            {
                // cerr << "!!!!!!!!!skip: " << item_i << endl;
                ++order_i;
                continue;
            }


            Vec2 next_pos(-1000000, -100000);

            Vec2 lower, upper;
            if (!range_to_get_item(cur_pos, items[item_i].region(), lower, upper))
            {
                // cur_posから直線に移動してitem_iを取れない場合
                int mask = intersect_holes_mask(cur_pos, items[item_i].pos(), holes);
                next_pos = next_dest_pos(cur_pos, cur_angle, items[item_i].pos(), mask);
                // cerr << cur_pos << " -> " << next_pos << endl;
                assert(!cur_pos.equals(next_pos));
            }
            else
            {
                // cur_posから角度調整すれば、直線に移動してitem_iが取れる
                // debug(cur_pos + lower);
                // debug(cur_pos + upper);
                assert(can_go_straight(cur_pos, cur_pos + lower, holes));
                assert(can_go_straight(cur_pos, cur_pos + upper, holes));
                assert(intersect(ori_items[item_i].region(), cur_pos, cur_pos + lower));
                assert(intersect(ori_items[item_i].region(), cur_pos, cur_pos + upper));


                // int item_j = item_i + 1;
                int order_j = order_i + 1;
                while (order_j < items.count())
                {
                    const int item_j = item_order[order_j];
                    if (intersect(items[item_j].region(), cur_pos))
                    {
                        ++order_j;
                        continue;
                    }
                    // debug(item_j);

                    // cerr << "one" << endl;
                    // printf("range: [%f, %f]\n", Math::RadToDeg(get_angle(lower)), Math::RadToDeg(get_angle(upper)));
                    Vec2 lower_j, upper_j;
                    if (!range_to_get_item(cur_pos, items[item_j].region(), lower_j, upper_j))
                        break;

                    // cerr << "two" << endl;
                    // printf("range: [%f, %f]\n", Math::RadToDeg(get_angle(lower)), Math::RadToDeg(get_angle(upper)));
                    Vec2 nlower, nupper;
                    if (!and_range(lower, upper, lower_j, upper_j, nlower, nupper))
                        break;
                    // cerr << "three" << endl;
                    // printf("range: [%f, %f]\n", Math::RadToDeg(get_angle(lower)), Math::RadToDeg(get_angle(upper)));

                    lower = nlower;
                    upper = nupper;

                    ++order_j;
                }

                // 進む方向を決める
                Vec2 move_vec;
                if (order_j == items.count())
                {
                    const Circle& last_item = items[item_order[items.count() - 1]].region();
                    lower = modify_out_lower(cur_pos, lower, last_item);
                    upper = modify_out_upper(cur_pos, upper, last_item);

                    Vec2 cur_vec(1, 0);
                    cur_vec.rotate(cur_angle);
                    if (lower.cross(cur_vec) < 0)
                    {
                        move_vec = lower;
                        for (;;)
                        {
                            float move_d = cur_pos.dist(intersect_point(items[item_i].region(), cur_pos, cur_pos + move_vec));
                            move_vec.normalize(move_d);

                            if (stage.field().isIn(cur_pos + move_vec))
                                break;
                            else
                                move_vec.rotate(Parameter::PLAYER_ANGLE_RATE);
                        }
                    }
                    else if (upper.cross(cur_vec) > 0)
                    {
                        move_vec = upper;
                        for (;;)
                        {
                            float move_d = cur_pos.dist(intersect_point(items[item_i].region(), cur_pos, cur_pos + move_vec));
                            move_vec.normalize(move_d);

                            if (stage.field().isIn(cur_pos + move_vec))
                                break;
                            else
                                move_vec.rotate(-Parameter::PLAYER_ANGLE_RATE);
                        }
                    }
                    else
                    {
                        move_vec = cur_vec;
                        float move_d = cur_pos.dist(intersect_point(items[item_i].region(), cur_pos, cur_pos + move_vec));
                        move_vec.normalize(move_d);
                    }
                }
                else
                {
                    const Vec2& ndest = next_dest[item_order[order_j - 1]];
                    const Vec2& ndest_vec = ndest - (order_j >= 2 ? items[item_order[order_j - 2]].pos() : cur_pos);
                    // move_vec = lower.cross(ndest_vec) < 0 ? lower : upper;
                    if (lower.cross(ndest_vec) < 0)
                    {
                        const int inf = 1000000;
                        int best_cost = inf;
                        Vec2 best_mvec;
                        Vec2 mvec = lower;
                        for (int i = 0; mvec.cross(upper) > 0 && (i < 30 || best_cost == inf); ++i)
                        {
                            float move_d = cur_pos.dist(intersect_point(items[item_i].region(),
                                                        cur_pos, cur_pos + mvec));
                            mvec.normalize(move_d);

                            int cost = calc_cost(cur_pos, cur_angle, cur_pos + mvec, ndest);
                            if (stage.field().isIn(cur_pos + mvec) && cost < best_cost)
                            {
                                best_cost = cost;
                                best_mvec = mvec;
                            }

                            mvec.rotate(Parameter::PLAYER_ANGLE_RATE / 2);
                        }
                        assert(best_cost < inf);
                        move_vec = best_mvec;
                    }
                    else
                    {
                        const int inf = 1000000;
                        int best_cost = inf;
                        Vec2 best_mvec;
                        Vec2 mvec = upper;
                        for (int i = 0; lower.cross(mvec) > 0 && (i < 30 || best_cost == inf); ++i)
                        {
                            float move_d = cur_pos.dist(intersect_point(items[item_i].region(),
                                                        cur_pos, cur_pos + mvec));
                            mvec.normalize(move_d);

                            int cost = calc_cost(cur_pos, cur_angle, cur_pos + mvec, ndest);
                            if (stage.field().isIn(cur_pos + mvec) && cost < best_cost)
                            {
                                best_cost = cost;
                                best_mvec = mvec;
                            }

                            mvec.rotate(-Parameter::PLAYER_ANGLE_RATE / 2);
                        }
                        assert(best_cost < inf);
                        move_vec = best_mvec;
                    }
                }
                // const float move_d = cur_pos.dist(items[item_i].pos());
                // cerr << endl;
                // debug(items[item_i].pos());
                // debug(items[item_i].radius());
                // debug(cur_pos);
                // debug(cur_pos + move_vec);
                next_pos = cur_pos + move_vec;


                ++order_i;
            }
            // printf("%s -> %s\n", to_s(cur_pos), to_s(next_pos));
            // assert(can_go_straight(cur_pos, next_pos, holes));

            // debug
            // if (true)
#ifdef LOCAL
            if (false)
            {
                // cerr << endl;
                // fprintf(stderr, "%s -> %s\n", to_s(cur_pos), to_s(next_pos));
                Circle cplayer(cur_pos, Parameter::PLAYER_RADIUS);
                Vec2 vec = next_pos - cplayer.pos();
                rep(i, items.count())
                {
                    if (debug_items[i].tryToRemove(cplayer, vec))
                    {
                        // cerr << "removed: " << i << endl;
                    }
                }

                if (order_i > 0)
                {
                    assert(debug_items[item_order[order_i - 1]].isRemoved());
                }
            }
#endif


            cur_angle = get_angle(cur_pos, next_pos);
            cur_pos = next_pos;
            turning_points[num_turning_points++] = cur_pos;
        }

        // cerr << "points: " << endl;
        // rep(i, num_turning_points)
        //     cerr << turning_points[i] << endl;

        return num_turning_points;
    }

    // helper methods for solve
private:
    Vec2 next_dest_pos(const Vec2& pos, const float angle, const Vec2& dest_item, int holes_mask)
    {
        assert(holes_mask);
        assert(!can_go_straight(pos, dest_item, holes));

        int n = 0;
        const Vec2* p[64];
        rep(i, holes.count())
        {
            if (holes_mask >> i & 1)
            {
                rep(j, 4)
                    p[n++] = &corners[i][j];
            }
        }
        const int start = n, goal = n + 1;
        p[start] = &pos;
        p[goal] = &dest_item;
        n += 2;

        const float inf = 1e7;
        int move_cost[64][64];
        float ang[64][64];
        rep(i, n) rep(j, n)
            move_cost[i][j] = inf;

        rep(i, n - 2) rep(j, n - 2)
        {
            if (i != j)
            {
                // FIXME: move_costの値がおかしい
                // move_cost[i][j] = corner_move_cost[i / 4][i % 4][j / 4][j % 4];
                move_cost[i][j] = need_move(p[i]->dist(*p[j])); // とりあえずその都度計算してる
                ang[i][j] = corner_angle[i / 4][i % 4][j / 4][j % 4];

                // if (move_cost[i][j] < 1000)
                //     cerr << *p[i] << " -> " << *p[j] << ": " << move_cost[i][j] << endl;
            } 
        }
        rep(i, n - 2) for (int j = start; j <= goal; ++j)
        {
            if (can_go_straight(*p[i], *p[j], holes) && !p[i]->equals(*p[j]))
            {
                move_cost[i][j] = move_cost[j][i] = need_move(p[i]->dist(*p[j]));
                ang[i][j] = get_angle(*p[i], *p[j]);
                ang[j][i] = get_angle(*p[j], *p[i]);
            }
        }

        const int unreach = 10000000;
        int dp[64][64]; // (cur, prev)
        int prev[64][64];
        rep(i, n) rep(j, n)
            dp[i][j] = unreach;
        PriorityQueue<uint, 64 * 64> q;
        rep(i, n)
        {
            if (move_cost[start][i] < inf)
            {
                dp[i][start] = move_cost[start][i] + need_rot(ang[start][i] - angle);
                q.push(pss(dp[i][start], pcc(i, start)));
            }
        }
        while (!q.empty())
        {
            uint tt = q.top();
            q.pop();

            int c = pss_first(tt);
            tt = pss_second(tt);
            int cur_pos = pcc_first(tt);
            int prev_pos = pcc_second(tt);

            if (cur_pos == goal)
            {
                int i = cur_pos, j = prev_pos;
                while (j != start)
                {
                    int ti = i, tj = j;
                    i = tj;
                    j = prev[ti][tj];
                }
                return *p[i];
            }

            rep(to, n)
            {
                if (move_cost[cur_pos][to] < inf)
                {
                    int nc = c + move_cost[cur_pos][to];
                    if (nc < dp[to][cur_pos]) // 角度計算の無駄を省くため2回に分けてる
                    {
                        const float cur_angle = ang[prev_pos][cur_pos];
                        const float next_angle = ang[cur_pos][to];
                        nc += need_rot(get_angle(cur_angle, next_angle));
                        if (nc < dp[to][cur_pos])
                        {
                            dp[to][cur_pos] = nc;
                            prev[to][cur_pos] = prev_pos;
                            q.push(pss(nc, pcc(to, cur_pos)));
                        }
                    }
                }
            }
        }
        abort();
        return Vec2(-1000000, -1000000);
    }

    // posからベクトル方向にdist進むと、穴に落っこちるときに使う
    // 穴に落っこちないような方向を求める
    // return == falseのときはlowerに変更なし
    bool modify_lower(const Vec2& pos, Vec2& lower, float dist, int holes_mask)
    {
        assert(holes_mask);

        bool ok = false;
        Vec2 mod_lower;
        rep(i, holes.count())
        {
            if (holes_mask >> i & 1)
            {
                rep(j, 4)
                {
                    const float eps = 0.01;
                    Vec2 vec = corners[i][j] - pos;
                    if (vec.length() > eps && lower.cross(vec) > 0)
                    {
                        vec *= dist / vec.length();
                        if (can_go_straight(pos, pos + vec, holes))
                        {
                            if (ok)
                            {
                                // より範囲が広くなるようならvecを使う
                                if (mod_lower.cross(vec) < 0)
                                    mod_lower = vec;
                            }
                            else
                                mod_lower = vec;

                            ok = true;
                        }
                    }
                }
            }
        }
        if (ok)
            lower = mod_lower;
        return ok;
    }
    bool modify_upper(const Vec2& pos, Vec2& upper, float dist, int holes_mask)
    {
        assert(holes_mask);

        bool ok = false;
        Vec2 mod_upper;
        rep(i, holes.count())
        {
            if (holes_mask >> i & 1)
            {
                rep(j, 4)
                {
                    const float eps = 0.01;
                    Vec2 vec = corners[i][j] - pos;
                    if (vec.length() > eps && upper.cross(vec) < 0)
                    {
                        vec *= dist / vec.length();
                        if (can_go_straight(pos, pos + vec, holes))
                        {
                            if (ok)
                            {
                                // より範囲が広くなるようならvecを使う
                                if (mod_upper.cross(vec) > 0)
                                    mod_upper = vec;
                            }
                            else
                                mod_upper = vec;

                            ok = true;
                        }
                    }
                }
            }
        }
        if (ok)
            upper = mod_upper;

        return ok;
    }

    // return: 有効な範囲があるか
    // lower, upperには接線ベクトルが入る
    bool range_to_get_item(const Vec2& cur_pos, const Circle& item, Vec2& lower, Vec2& upper)
    {
        tangent_vecs(cur_pos, item, lower, upper);
        const float dist_to_item = lower.length();
        const int lower_mask = intersect_holes_mask(cur_pos, cur_pos + lower, holes);
        const int upper_mask = intersect_holes_mask(cur_pos, cur_pos + upper, holes);

        if (lower_mask)
            if (!modify_lower(cur_pos, lower, dist_to_item, lower_mask))
                return false;
        if (upper_mask)
            if (!modify_upper(cur_pos, upper, dist_to_item, upper_mask))
                return false;
        if (lower.cross(upper) < 0)
            return false;
        return true;
    }

    // 円が場外にはみ出している時、stage外に出ないlowerを求める
    Vec2 modify_out_lower(const Vec2& cur_pos, const Vec2& cur_lower, const Circle& c)
    {
        const Rectangle& field = expand_margin(stage.field(), -1e-4); // 誤差死回避のために少し小さくする

        if (field.isIn(intersect_point(c, cur_pos, cur_pos + cur_lower)))
            return cur_lower;

        Vec2 best_lower = c.pos() - cur_pos;
        Vec2 corner[] = { field.pointA(), field.pointB(), field.pointC(), field.pointD() };
        rep(i, 4)
        {
            const Vec2& a = corner[i], b = corner[(i + 1) % 4];
            if (intersect(c, a, b))
            {
                Vec2 ip1, ip2;
                intersect_points(c, a, b, ip1, ip2);
                Vec2 vec1 = ip1 - cur_pos, vec2 = ip2 - cur_pos;
                if (vec1.cross(cur_lower) < 0 && best_lower.cross(vec1) < 0)
                    best_lower = vec1;
                if (vec2.cross(cur_lower) < 0 && best_lower.cross(vec2) < 0)
                    best_lower = vec2;
            }
        }
        return best_lower;
    }
    Vec2 modify_out_upper(const Vec2& cur_pos, const Vec2& cur_upper, const Circle& c)
    {
        const Rectangle& field = expand_margin(stage.field(), -1e-4); // 誤差死回避のために少し小さくする

        if (field.isIn(intersect_point(c, cur_pos, cur_pos + cur_upper)))
            return cur_upper;

        Vec2 best_upper = c.pos() - cur_pos;
        Vec2 corner[] = { field.pointA(), field.pointB(), field.pointC(), field.pointD() };
        rep(i, 4)
        {
            const Vec2& a = corner[i], b = corner[(i + 1) % 4];
            if (intersect(c, a, b))
            {
                Vec2 ip1, ip2;
                intersect_points(c, a, b, ip1, ip2);
                Vec2 vec1 = ip1 - cur_pos, vec2 = ip2 - cur_pos;
                if (vec1.cross(cur_upper) > 0 && best_upper.cross(vec1) > 0)
                    best_upper = vec1;
                if (vec2.cross(cur_upper) > 0 && best_upper.cross(vec2) > 0)
                    best_upper = vec2;
            }
        }
        return best_upper;
    }
private:

    const Stage& stage;
    ItemCollection items;
    ItemCollection ori_items;
    HoleCollection holes;

    Vec2 corners[16][4];
    float corner_move_cost[16][4][16][4];
    float corner_angle[16][4][16][4];
};
void make_actions(const Stage& stage, const Vec2* turning_points, int num_turning_points, Answer& ans)
{
    ans.reset();

    Player player = stage.player();
    // rep(i, num_turning_points)
    for (int i = 0; i < num_turning_points; )
    {
        const Vec2& dest = turning_points[i];
        if (player.pos().equals(dest))
        {
            ++i;
            continue;
        }

        // const Vec2 vec_to_item = dest - player.pos();
        // const Vec2 player_dir = player.vec();
        const float dest_angle = get_angle(player.pos(), dest);
        const float cur_angle = player.arg();
        const float diff_angle = normalize_angle(dest_angle - cur_angle);
        // if (abs(player_dir.cos(vec_to_item) - 1) > 1e-3)
        if (abs(diff_angle) > 1e-4)
        {
            // float rot = Math::LimitAbs(player_dir.rotSign(vec_to_item), Parameter::PLAYER_ANGLE_RATE);
            float rot = Math::LimitAbs(diff_angle, Parameter::PLAYER_ANGLE_RATE);
            hpc::Action act(ActionType_Rotate, rot);
            ans.add(act);
            player.act(act);
        }
        else
        {
            // debug(to_s(player.pos()));
            // assert(player.pos().x < 55);

            const float dist = player.pos().dist(dest);
            const float move_d = min(dist, Parameter::PLAYER_SPEED);
            hpc::Action act(ActionType_Move, move_d);
            ans.add(act);
            player.act(act);
        }
    }
}

int calc_cost_call = 0;
class OrderSolver
{
public:
    OrderSolver(const Stage& _stage)
        : stage(_stage), items(stage.items()), n(items.count()), player(_stage.player())
    {
        rep(i, n)
            order[i] = i;


        rep(i, n)
            pos[i] = items[i].pos();
        pos[n] = player.pos();

        rep(j, n + 1) rep(i, j + 1)
            mcost[i][j] = mcost[j][i] = need_move(pos[i].dist(pos[j]));

        rep(i, n + 1) rep(j, n + 1)
            ang[i][j] = get_angle(pos[j] - pos[i]);
    }

    void get_order(int* res) const
    {
        rep(i, n)
            res[i] = order[i];
    }

    void solve()
    {
        greedy_near();
        improve_simply();
    }

    int calc_move_cost() const
    {
        int move_cost = 0;
        move_cost += mcost[n][order[0]];
        rep(i, n - 1)
            move_cost += mcost[order[i]][order[i + 1]];
        return move_cost;
    }
    int calc_rot_cost() const
    {
        int rot_cost = 0;
        rot_cost += need_rot(get_angle(ang[n][order[0]], player.arg()));
        float angle = ang[n][order[0]];
        rep(i, n - 1)
        {
            const float next_angle = ang[order[i]][order[i + 1]];
            rot_cost += need_rot(get_angle(angle, next_angle));
            angle = next_angle;
        }
        return rot_cost;
    }
    int calc_cost() const
    {
        // move_costを2倍するとなぜか伸びる
        return 2 * calc_move_cost() + calc_rot_cost();
    }
private:
    void greedy_near()
    {
        bool used[128] = {};

        Vec2 p = stage.player().pos();
        // float a = stage.player().arg();
        rep(oi, n)
        {
            int k = -1;
            float min_cost = 1e7;
            rep(i, n)
            {
                if (!used[i])
                {
                    float d = p.dist(items[i].pos());
                    // float rot = get_angle(p, a, items[i].pos());
                    // float c = d + rot;
                    float c = d;
                    if (c < min_cost)
                    {
                        min_cost = c;
                        k = i;
                    }
                }
            }
            assert(k != -1);

            order[oi] = k;
            p = items[k].pos();
            // a = get_angle(items[k].pos() - p);
            used[k] = true;
        }
    }

    void improve_simply()
    {
        if (n < 3)
            return;

        static Random rand;
        int cur = calc_cost();
        int ori[128];
        rep(i, n)
            ori[i] = order[i];
        // for (int loop = 0; loop < n * n * 10; ++loop)
        for (;;)
        {
            bool updated = false;
            rep(b, n) rep(a, b - 1)
            {
                // int a = rand.next_int(n - 2);
                // int b = a + rand.next_int(n - a - 1);

                swap(order[a + 1], order[b]);
                reverse(order + a + 2, order + b);

                int next = calc_cost();
                if (next < cur)
                {
                    updated = true;
                    cur = next;
                    rep(i, n)
                        ori[i] = order[i];
                }
                else
                {
                    rep(i, n)
                        order[i] = ori[i];
                }
            }
            if (!updated)
                break;
        }
    }


private:
    const Stage& stage;
    const ItemCollection& items;
    const int n;
    const solver::Player player;

    int order[128];

    Vec2 pos[128];
    int mcost[128][128];
    float ang[128][128];
};

#ifdef LOCAL
int expected_move, expected_rot;
#endif
void solve_order(const Stage& stage, int* order)
{
    OrderSolver s(stage);
    s.solve();
    s.get_order(order);

#ifdef LOCAL
    expected_move = s.calc_move_cost();
    expected_rot = s.calc_rot_cost();
#endif
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
                if (mi == num_mark - 1 && oi < n &&
                    // IsHit(item.region(), player.region(), cur + cur_vec) &&
                    titem.tryToRemove(player.region(), 100 * cur_vec) &&
                    (oi == n - 1 || cur_vec.cross(dest - cur) * cur_vec.cross(items[oi + 1].pos() - cur) < 0) &&
                    stage.field().isIn(cur + d * cur_vec) &&
                    can_go_straight(cur, cur + d * cur_vec, holes)
                   )
                {
                    assert(dest.equals(item.pos()));
                    force_move = true;
                    // cerr << "force" << endl;
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
            // cerr << "ouou" << endl;
        }
        // cerr << "loopend" << endl;
    }
}


Answer answer;
int ans_i;
void solve(const Stage& stage)
{
    int order[128];
    solve_order(stage, order);
    // solve_actions(stage, order, answer);


    ActionSolver sol(stage);
    Vec2 pos[1024];
    // cerr << "in solve" << endl;
    int n = sol.solve(order, pos);
    // cerr << "out solve" << endl;
    // debug(n);
    // rep(i, n)
    //     fprintf(stderr, "(%3.3f, %3.3f)\n", pos[i].x, pos[i].y);

    solver::make_actions(stage, pos, n, answer);

    ans_i = 0;
}

hpc::Action get_next_action()
{
    const Action& t = answer.action(ans_i++);
    return hpc::Action(t.type(), t.value());
}

double get_score(const Stage& stage)
{
    int n = stage.items().count();
    double score = (double)n * n / (answer.num_move() + answer.num_rot()) * 10000.0;
    return score;
}
#ifdef LOCAL
void print_score(const Stage& stage)
{
    fprintf(stderr, "%4d %4d %10.3f ", answer.num_move(), answer.num_rot(), get_score(stage));
}
void print_info()
{
    fprintf(stderr, "%4d (%+4d), %4d (%+4d)\n", answer.num_move(), answer.num_move() - expected_move,answer.num_rot(), answer.num_rot() - expected_rot);
}
#endif

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
        // if (solver::stage_no > 35)
        //     exit(1);
        // cerr << endl;
        // cerr << string(60, '=') << endl;
        // debug(solver::stage_no);

        solver::solve(aStage);
        // printf("%d, %d\n", solver::answer.num_move(), solver::answer.num_rot());
        // debug(solver::na);
        // debug(solver::calc_cost_call);
        // static int move = 0, rot = 0;
        // printf("%4d %4d (%.3f) %8d %8d\n",
        //         solver::answer.num_move(), solver::answer.num_rot(),
        //         float(solver::answer.num_rot()) / solver::answer.num_move(),
        //         move, rot);

        // solver::print_score(aStage);
        // solver::print_info();

        // static float sum = 0;
        // sum += solver::get_score(aStage);
        // printf("sum: %f\n", sum);
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
        if (solver::ans_i == solver::answer.size())
            return Action(ActionType_Rotate, 1);
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

        // Player player = aStage.player();
        // printf("%s, %f\n", solver::to_s(aStage.player().pos()).c_str(), Math::RadToDeg(aStage.player().arg()));

        // ItemCollection items = aStage.items();
        // vector<int> remain;
        // rep(i, items.count())
        //     if (!items[i].isRemoved())
        //         remain.push_back(i);
        // debug(remain);

        return solver::get_next_action();
    }
}

//----------------------------------------------------------
// EOF
