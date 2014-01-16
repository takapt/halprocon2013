#include "HPCAnswerInclude.hpp"

#undef LOCAL
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
#define debug(a)
#endif

#define RUNTIME_DEBUG


#define rep(i, n) for (int i = 0; i < (int)(n); ++i)

typedef unsigned int uint;
typedef unsigned short ushort;
typedef unsigned char uchar;

namespace solver
{
using namespace hpc;

#ifdef LOCAL
using namespace std;

float get_angle(const Vec2& vec);
float to_deg(const Vec2& vec)
{
    return Math::RadToDeg(get_angle(vec));
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

//////////////////////////////////
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
// 円cと直線abとの交点
// return: ip1, ip2
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
// 円cと線分(from, to)
bool intersect(const Circle& c, const Vec2& from, const Vec2& to)
{
    // 動いていなければ静止している円の判定
    if (from == to)
        return intersect(c, from);

    // 実質的に半径 aC0.radius + aC1.radius の円と線分 [aC1.pos, aP1.pos] の衝突判定と
    // 同等になる。
    const Vec2 seg = to - from;
    const Vec2 from_to_c = c.pos() - from;
    const float dist = Math::Abs(seg.cross(from_to_c)) / seg.length();
    // 距離が c.radius より遠ければ衝突しない
    if(dist > c.radius()) {
        return false;
    }
    // 線分の延長線上で交差していないか調べる。
    const Vec2 to_to_c = c.pos() - to;
    // それぞれの点が円の反対方向にあれば衝突
    if (from_to_c.dot(seg) * to_to_c.dot(seg) <= 0)
        return true;
    // 半径が大きければ衝突
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


// hpc::Actionに引数なしのコンストラクタがなかったので
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

// [-PI, PI]?
inline bool valid_angle(float angle)
{
    return -Math::PI <= angle && angle <= Math::PI;
}
// angleを[-PI, PI]の範囲にする
inline float normalize_angle(float angle)
{
    assert(!isnan(angle));
    assert(!isinf(angle));

    while (angle > Math::PI)
        angle -= 2 * Math::PI;
    while (angle < -Math::PI)
        angle += 2 * Math::PI;

    assert(valid_angle(angle));
    return angle;
}
// x軸と(pos.x, pos.y)のなす角度
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
// angleをdest_angleにするために必要な回転角度
// 関数名がわかりにくい...orz
inline float get_angle(float angle, float dest_angle)
{
    return normalize_angle(dest_angle - angle);
}
// プレイヤーが座標pos、角度angleのときに、dest方向へ向くために必要な回転角度
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
        return (int)(dist) + 1;
}
inline int need_rot(float angle)
{
    angle = normalize_angle(angle);
    assert(valid_angle(angle));

    if (angle < 0)
        angle = -angle;
    const float eps = 1e-7;
    if (angle < eps)
        return 0;
    else
        return (int)(angle * 10) + 1;
}

// cur -> pos[0] -> pos[1] -> ... -> pos[n - 1]と移動するときに必要なターン数
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
// cur -> a -> bと行動するときに必要なターン数
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

    double score(const Stage& stage) const
    {
        const int items = stage.items().count();
        return items * items / double(num_move() + num_rot()) * 10000;
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


// 線分abと線分cdが交差するか
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
// 線分abと長方形rectが交差するか
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

// 高速化アルゴリズム適応前に使っていた愚直な方法
bool _can_go_straight(const Vec2& a, const Vec2& b, const HoleCollection& holes)
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

// pointから引いた、circleとの接線
// return: pointからの接線ベクトルlower, upper
// post condition: lower.cross(upper) > 0
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

    assert(lower.cross(upper) > 0);
}

// 高速化アルゴリズム適応前に使っていた愚直な方法
// return: hole iと交わる時、result |= 1 << i
int _intersect_holes_mask(const Vec2& from, const Vec2& to, const HoleCollection& holes)
{
    int mask = 0;
    rep(i, holes.count())
    {
        if (intersect(from, to, holes[i]))
            mask |= 1 << i;
    }
    return mask;
}

// 角度の範囲[lower_a, upper_a]と[lower_b, upper_b]の共通部分を返す
// return:
//  返り値有効な範囲があるか
//  共通部分はlower, upperに入る
bool and_range(const Vec2& lower_a, const Vec2& upper_a, const Vec2& lower_b, const Vec2& upper_b, Vec2& lower, Vec2& upper)
{
    if (lower_a.cross(upper_b) < 0 || upper_a.cross(lower_b) > 0)
        return false;

    lower = lower_a.cross(lower_b) > 0 ? lower_b : lower_a;
    upper = upper_a.cross(upper_b) < 0 ? upper_b : upper_a;

    return lower.cross(upper) > 0;
}

// 穴に関する処理を高速化するためのデータ構造
class HoleManager
{
public:
    HoleManager() : n(-1000), w(-1000), h(-100){}

    HoleManager(float _w, float _h, const HoleCollection& _holes)
        : n(_holes.count()), w(_w + 1e-4), h(_h + 1e-4)
    {
        rep(i, n)
        {
            holes[i] = _holes[i];
            corners[i][0] = holes[i].pointA();
            corners[i][1] = holes[i].pointB();
            corners[i][2] = holes[i].pointC();
            corners[i][3] = holes[i].pointD();
        }

        rep(y, h + 1) rep(x, w + 1)
            holes_mask[y][x] = 0;
        rep(i, w + 1) rep(j, w + 1)
            hor_mask[i][j] = 0;
        rep(i, h + 1) rep(j, h + 1)
            ver_mask[i][j] = 0;

        rep(hi, n)
        {
            const uint bit = 1 << hi;

            const int left = holes[hi].left;
            const int right = holes[hi].right + 1e-4;
            const int bottom = holes[hi].bottom;
            const int top = holes[hi].top + 1e-4;
            assert(0 <= left);
            assert(right < w);
            assert(0 <= bottom);
            assert(top < h);

            for (int y = bottom; y <= top; ++y)
                for (int x = left; x <= right; ++x)
                    holes_mask[y][x] |= bit;

            for (int i = 0; i <= right; ++i)
                for (int j = left; j < w; ++j)
                    hor_mask[i][j] |= bit;

            for (int i = 0; i <= top; ++i)
                for (int j = bottom; j < h; ++j)
                    ver_mask[i][j] |= bit;
        }
    }

    HoleManager(float _w, float _h, Rectangle* _holes, int num)
        : n(num), w(_w + 1e-4), h(_h + 1e-4)
    {
        rep(i, n)
        {
            holes[i] = _holes[i];
            corners[i][0] = holes[i].pointA();
            corners[i][1] = holes[i].pointB();
            corners[i][2] = holes[i].pointC();
            corners[i][3] = holes[i].pointD();
        }

        rep(y, h + 1) rep(x, w + 1)
            holes_mask[y][x] = 0;
        rep(i, w + 1) rep(j, w + 1)
            hor_mask[i][j] = 0;
        rep(i, h + 1) rep(j, h + 1)
            ver_mask[i][j] = 0;

        rep(hi, n)
        {
            const uint bit = 1 << hi;

            const int left = holes[hi].left;
            const int right = holes[hi].right + 1e-4;
            const int bottom = holes[hi].bottom;
            const int top = holes[hi].top + 1e-4;
            assert(0 <= left);
            assert(right < w);
            assert(0 <= bottom);
            assert(top < h);

            for (int y = bottom; y <= top; ++y)
                for (int x = left; x <= right; ++x)
                    holes_mask[y][x] |= bit;

            for (int i = 0; i <= right; ++i)
                for (int j = left; j < w; ++j)
                    hor_mask[i][j] |= bit;

            for (int i = 0; i <= top; ++i)
                for (int j = bottom; j < h; ++j)
                    ver_mask[i][j] |= bit;
        }
    }


    // 線分abが交差している可能性のある穴はbitを立てる
    // bitが立っている場合は、交差している or していない
    // bitが立っていない場合は、交差していない
    uint possible_mask(const Vec2& a, const Vec2& b) const
    {
        float _left = a.x, _right = b.x;
        float _bottom = a.y, _top = b.y;
        if (_left > _right)
            swap(_left, _right);
        if (_bottom > _top)
            swap(_bottom, _top);

        int left = max<int>(0, _left), right = min<int>(w - 1, _right + 1e-4);
        int bottom = max<int>(0, _bottom), top = min<int>(h - 1, _top + 1e-4);
        assert(0 <= left);
        assert(right < w);
        assert(0 <= bottom);
        assert(top < h);

        return hor_mask[left][right] & ver_mask[bottom][top];
    }

    uint intersect_holes_mask(const Vec2& a, const Vec2& b) const
    {
        // 可能性のある穴とだけ交差判定をすることで高速化
        uint res = 0;
        int possi = possible_mask(a, b);
        rep(i, n)
        {
            if (possi >> i & 1)
            {
                if (intersect(a, b, i))
                    res |= 1 << i;
            }
        }
        return res;
    }

    bool intersect(const Vec2& a, const Vec2& b, int hole_i) const
    {
        assert(0 <= hole_i && hole_i < n);
        rep(i, 4)
        {
            if (solver::intersect(a, b, corners[hole_i][i], corners[hole_i][(i + 1) % 4]))
                return true;
        }
        return false;
    }

    bool can_go_straight(const Vec2& a, const Vec2& b, uint mask) const
    {
        rep(i, n)
        {
            if (mask >> i & 1)
            {
                if (intersect(a, b, i))
                    return false;
            }
        }
        return true;
    }
    bool can_go_straight(const Vec2& a, const Vec2& b) const
    {
        return can_go_straight(a, b, possible_mask(a, b));
    }

    int size() const { return n; }
private:
    int n;
    int w, h;
    Rectangle holes[16];
    Vec2 corners[16][4];

    uint holes_mask[64][64];
    uint hor_mask[64][64];
    uint ver_mask[64][64];
};
HoleManager hole_manager;
bool can_go_straight(const Vec2& a, const Vec2& b, const HoleCollection& holes)
{
    return hole_manager.can_go_straight(a, b);
}
int intersect_holes_mask(const Vec2& from, const Vec2& to, const HoleCollection& holes)
{
    const int res = hole_manager.intersect_holes_mask(from, to);
    return res;
}

// アイテムに関する処理を高速化するためのデータ構造
class ItemManager
{
public:
    ItemManager() {}

    ItemManager(float _w, float _h, ItemCollection _items)
        : w(_w + 1e-4), h(_h + 1e-4), items(_items), m(0)
    {
        int ii = 0;
        while (ii < items.count())
        {
            Rectangle rects[SIZE];
            int nr = 0;
            for (int i = 0; i < SIZE && ii < items.count(); ++i, ++ii)
            {
                const Circle& c = items[ii].region();
                const float r = c.radius();
                const float x = c.pos().x, y = c.pos().y;
                const float left = max(0.f, x - r);
                const float right = min(w, x + r);
                const float bottom = max(0.f, y - r);
                const float top = min(h, y + r);

                rects[nr++] = Rectangle(left, right, bottom, top);
            }

            hole_m[m++] = HoleManager(w + 5, h + 5, rects, nr);
        }
    }

    int manager_index(int i) const
    {
        return i / SIZE;
    }
    int index_in_manager(int i) const
    {
        return i % SIZE;
    }

    bool intersect(const Vec2& a, const Vec2& b, int item_i) const
    {
        const HoleManager& hm = hole_m[manager_index(item_i)];
        int i = index_in_manager(item_i);

        int possi_mask = hm.possible_mask(a, b);
        if (possi_mask >> i & 1)
            return solver::intersect(items[item_i].region(), a, b);
        else
            return false;
    }

private:
    float w, h;
    ItemCollection items;
    int m;
    HoleManager hole_m[4];
    static const int SIZE = 32;
};
ItemManager item_manager;


// アイテムの取得順を与えて、アクションを決定する
class ActionSolver
{
public:
    ActionSolver(const Stage& _stage)
        : stage(_stage), holes(_stage.holes())
    {
        const float w = stage.field().width();
        const float h = stage.field().height();
        hole_manager = HoleManager(w, h, stage.holes());

        const float hole_margin = 0.001;
        rep(i, holes.count())
        {
            Rectangle rect = expand_margin(holes[i], hole_margin);
            corners[i][0] = rect.pointA();
            corners[i][1] = rect.pointB();
            corners[i][2] = rect.pointC();
            corners[i][3] = rect.pointD();
        }

        const ItemCollection& _items = stage.items();
        rep(i, _items.count())
        {
            static const float eps = 1e-2;
            items.setupAddItem(_items[i].pos(), _items[i].radius() + Parameter::PLAYER_RADIUS - eps);
            ori_items.setupAddItem(_items[i].pos(), _items[i].radius() + Parameter::PLAYER_RADIUS);
        }
        item_manager = ItemManager(w, h, ori_items);
    }

    // return:
    //  返り値はturning_pointsの要素数
    //  turning_pointsの各要素はプレイヤーが回転する座標
    int solve(int* item_order, Vec2* turning_points)
    {
#ifdef LOCAL
        ItemCollection debug_items = stage.items();
#endif

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

        // got[i] -> アイテムiは取得済みか
        bool got[128];
        rep(i, items.count())
            got[i] = false;

        int num_turning_points = 0;
        turning_points[num_turning_points++] = cur_pos;
        for (int order_i = 0; order_i < items.count(); )
        {
            // 今回取るべきアイテムはitem_i
            const int item_i = item_order[order_i];

            if (got[item_i])
            {
                ++order_i;
                continue;
            }

            if (intersect(ori_items[item_i].region(), cur_pos))
            {
                ++order_i;
                continue;
            }


            Vec2 next_pos(-1000000, -100000); // 次に向かう座標

            // 次に向かう方向は[lower, upper]の範囲内でないといけない
            Vec2 lower, upper;
            if (!range_to_get_item(cur_pos, items[item_i].region(), lower, upper))
            {
                // cur_posから直線に移動してitem_iを取れない場合
                int mask = intersect_holes_mask(cur_pos, items[item_i].pos(), holes);
                tangent_vecs(cur_pos, items[item_i].region(), lower, upper);
                mask |= intersect_holes_mask(cur_pos, cur_pos + lower, holes);
                mask |= intersect_holes_mask(cur_pos, cur_pos + upper, holes);

                // next_posは穴の角になっている
                next_pos = next_dest_pos(cur_pos, cur_angle, items[item_i].pos(), mask);

                assert(!cur_pos.equals(next_pos));
            }
            else
            {
                // cur_posから角度調整すれば、直線に移動してitem_iが取れる
                assert(intersect(ori_items[item_i].region(), cur_pos, cur_pos + lower.getNormalized(100)));
                assert(intersect(ori_items[item_i].region(), cur_pos, cur_pos + upper.getNormalized(100)));


                // どのアイテムまでcur_posから直線に移動して取得可能か
                // 最終的に求まったorder_jを用いると、
                // order_i, order_i + 1, ..., order_j - 1のアイテムが直線移動で取得可能
                // ただし、直線移動する向きは[lower, upper]の範囲でないといけない
                int order_j = order_i + 1;
                while (order_j < items.count())
                {
                    const int item_j = item_order[order_j];
                    if (got[item_j] || intersect(items[item_j].region(), cur_pos))
                    {
                        ++order_j;
                        continue;
                    }

                    Vec2 lower_j, upper_j;
                    if (!range_to_get_item(cur_pos, items[item_j].region(), lower_j, upper_j))
                        break;

                    Vec2 nlower, nupper;
                    if (!and_range(lower, upper, lower_j, upper_j, nlower, nupper))
                        break;

                    lower = nlower;
                    upper = nupper;
                    assert(intersect(ori_items[item_i].region(), cur_pos, cur_pos + lower.getNormalized(100)));
                    assert(intersect(ori_items[item_i].region(), cur_pos, cur_pos + upper.getNormalized(100)));

                    ++order_j;
                }


                // 進む方向を[lower, upper]の範囲内で決める
                Vec2 move_vec;
                if (order_j == items.count())
                {
                    // 最後のアイテム
                    // 最後なので無駄回転する必要はない
#ifdef LOCAL
                    const Circle& last_ori_item = ori_items[item_order[items.count() - 1]].region();
                    if (!got[item_order[order_j - 1]] && !intersect(last_ori_item, cur_pos))
                    {
                        assert(intersect(last_ori_item, cur_pos, cur_pos + lower.getNormalized(100)));
                        assert(intersect(last_ori_item, cur_pos, cur_pos + upper.getNormalized(100)));
                    }
#endif

                    assert(intersect(ori_items[item_i].region(), cur_pos, cur_pos + lower.getNormalized(100)));
                    assert(intersect(ori_items[item_i].region(), cur_pos, cur_pos + upper.getNormalized(100)));

                    assert(lower.cross(upper) > 0);
#ifdef LOCAL
                    for (int i = order_i; i < order_j; ++i)
                    {
                        if (!got[item_order[i]] && !intersect(ori_items[item_order[i]].region(), cur_pos))
                        {
                            assert(intersect(ori_items[item_order[i]].region(), cur_pos, cur_pos + lower.getNormalized(100)));
                            assert(intersect(ori_items[item_order[i]].region(), cur_pos, cur_pos + upper.getNormalized(100)));
                        }
                    }
#endif

                    Vec2 cur_vec(1, 0);
                    cur_vec.rotate(cur_angle);
                    if (lower.cross(cur_vec) < 0)
                    {
                        move_vec = lower;
                        while (move_vec.cross(upper) >= 0)
                        {
                            float move_d = cur_pos.dist(intersect_point(items[item_i].region(), cur_pos, cur_pos + move_vec));
                            move_d = max(move_d, 1e-4f);
                            move_vec.normalize(move_d);

                            if (stage.field().isIn(cur_pos + move_vec))
                                break;
                            else
                                move_vec.rotate(Parameter::PLAYER_ANGLE_RATE / 10);
                        }
                    }
                    else if (upper.cross(cur_vec) > 0)
                    {
                        move_vec = upper;
                        while (lower.cross(move_vec) >= 0)
                        {
                            float move_d = cur_pos.dist(intersect_point(items[item_i].region(), cur_pos, cur_pos + move_vec));
                            move_d = max(move_d, 1e-4f);
                            move_vec.normalize(move_d);

                            if (stage.field().isIn(cur_pos + move_vec))
                                break;
                            else
                                move_vec.rotate(-Parameter::PLAYER_ANGLE_RATE / 10);
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
                    if (lower.cross(ndest_vec) < 0)
                    {
                        const int inf = 1000000;
                        int best_cost = inf;
                        Vec2 best_mvec;
                        Vec2 mvec = lower;
                        for (int i = 0; mvec.cross(upper) > 0 && (i < 20 || best_cost == inf); ++i)
                        {
                            float move_d = cur_pos.dist(intersect_point(items[item_i].region(),
                                                        cur_pos, cur_pos + mvec));

                            if (move_d < 1e-6)
                                break;
                            mvec.normalize(move_d);

                            int cost = calc_cost(cur_pos, cur_angle, cur_pos + mvec, ndest);
                            if (stage.field().isIn(cur_pos + mvec) && cost < best_cost)
                            {
                                best_cost = cost;
                                best_mvec = mvec;
                            }

                            mvec.rotate(Parameter::PLAYER_ANGLE_RATE / 2);
                        }
                        if (best_cost == inf)
                        {
                            best_cost = -1;
                            best_mvec = upper;
                        }
                        assert(best_cost < inf);
                        move_vec = best_mvec;
                    }
                    else
                    {
                        const int inf = 1000000;
                        Vec2 mvec = upper;
                        int best_cost = inf;
                        Vec2 best_mvec;
                        for (int i = 0; lower.cross(mvec) > 0 && (i < 20 || best_cost == inf); ++i)
                        {
                            float move_d = cur_pos.dist(intersect_point(items[item_i].region(),
                                                        cur_pos, cur_pos + mvec));
                            if (move_d < 1e-6)
                                break;
                            mvec.normalize(move_d);

                            int cost = calc_cost(cur_pos, cur_angle, cur_pos + mvec, ndest);
                            if (stage.field().isIn(cur_pos + mvec) && cost < best_cost)
                            {
                                best_cost = cost;
                                best_mvec = mvec;
                            }

                            mvec.rotate(-Parameter::PLAYER_ANGLE_RATE / 2);
                        }
                        if (best_cost == inf)
                        {
                            best_cost = -1;
                            best_mvec = lower;
                        }
                        assert(best_cost < inf);
                        move_vec = best_mvec;
                    }
                }
                next_pos = cur_pos + move_vec;

                ++order_i;
            }


            rep(i, items.count())
            {
#ifdef LOCAL
                bool fi = intersect(ori_items[i].region(), cur_pos, next_pos);
                bool fm = item_manager.intersect(cur_pos, next_pos, i);
                assert(fi == fm);
#endif
                if (!got[i] && item_manager.intersect(cur_pos, next_pos, i))
                {
                    got[i] = true;
                }
            }

            cur_angle = get_angle(cur_pos, next_pos);
            cur_pos = next_pos;
            turning_points[num_turning_points++] = cur_pos;
        }

        return num_turning_points;
    }

private:
    Vec2 next_dest_pos(const Vec2& pos, const float angle, const Vec2& dest_item, int holes_mask)
    {
        Vec2 res[128];
        list_turning_points(pos, angle, dest_item, holes_mask, res);
        return res[0];
    }
    // pos -> dest_itemに直線で行けないときに使う
    // posからdest_itemに行くまでに通過すべき点を列挙
    // 返り値num_pはresの要素数
    // pos -> (res[0] -> res[1] -> .. -> res[num_p - 1]) -> dest_item
    // 上記の()内が結果
    int list_turning_points(const Vec2& pos, const float angle, const Vec2& dest_item, int holes_mask, Vec2* res)
    {
        // holes_maskのbitが立っている穴の角を頂点として、dijkstraしている
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
                move_cost[i][j] = need_move(p[i]->dist(*p[j]));
                ang[i][j] = get_angle(*p[i], *p[j]);
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
                int num_p = 0;
                int i = cur_pos, j = prev_pos;
                while (j != start)
                {
                    int ti = i, tj = j;
                    i = tj;
                    j = prev[ti][tj];
                }
                res[num_p++] = *p[i];
                reverse(res, res + num_p);
                return num_p;
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
        debug(holes_mask);
        debug(pos);
        debug(dest_item);
        abort();
        return -1;
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
    // lower, upperに結果が入る
    // cur_posから[lower, upper]の範囲内の向きで直線に移動するとitemが取れる
    bool range_to_get_item(const Vec2& cur_pos, const Circle& item, Vec2& lower, Vec2& upper)
    {
        tangent_vecs(cur_pos, item, lower, upper);
        const float dist_to_item = lower.length();
        const int lower_mask = intersect_holes_mask(cur_pos, cur_pos + lower, holes);
        const int upper_mask = intersect_holes_mask(cur_pos, cur_pos + upper, holes);

        if (lower_mask)
        {
            if (!modify_lower(cur_pos, lower, dist_to_item, lower_mask))
                return false;
            if (lower.cross(upper) < 0)
                return false;
        }
        if (upper_mask)
        {
            if (!modify_upper(cur_pos, upper, dist_to_item, upper_mask))
                return false;
            if (lower.cross(upper) < 0)
                return false;
        }

        if (!intersect(item, cur_pos))
        {
            lower = modify_out_lower(cur_pos, lower, item);
            upper = modify_out_upper(cur_pos, upper, item);
        }
        return lower.cross(upper) > 0;
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

    // 半径がr + PLAYER_RADIUS - eps
    // -epsしているのは、アイテムとの交差がギリギリのときに、実際は取れてないというケースを防ぐため
    ItemCollection items;

    // 半径がr + PLAYER_RADIUS
    ItemCollection ori_items;

    HoleCollection holes;

    Vec2 corners[16][4]; // 穴の角
};

// turning_points[0] -> turning_points[1] -> ... -> turning_points[num_turning_points - 1]と移動する時に
// 必要なアクションを求める
void make_actions(const Stage& stage, const Vec2* turning_points, int num_turning_points, Answer& ans)
{
    ans.reset();

    Player player = stage.player();
    for (int i = 0; i < num_turning_points; )
    {
        const Vec2& dest = turning_points[i];
        if (player.pos().equals(dest))
        {
            ++i;
            continue;
        }

        const float dest_angle = get_angle(player.pos(), dest);
        const float cur_angle = player.arg();
        const float diff_angle = normalize_angle(dest_angle - cur_angle);
        if (abs(diff_angle) > 1e-4)
        {
            float rot = Math::LimitAbs(diff_angle, Parameter::PLAYER_ANGLE_RATE);
            hpc::Action act(ActionType_Rotate, rot);
            ans.add(act);
            player.act(act);
        }
        else
        {
            const float dist = player.pos().dist(dest);
            const float move_d = min(dist, Parameter::PLAYER_SPEED);
            hpc::Action act(ActionType_Move, move_d);
            ans.add(act);
            player.act(act);
        }
    }
}

int calc_cost(const Stage& stage, const Vec2* turning_points, int n)
{
    int total_cost = 0;
    Vec2 pos = stage.player().pos();
    float angle = stage.player().arg();
    rep(i, n)
    {
        float dest_angle = get_angle(pos, turning_points[i]);

        total_cost += need_move(pos.dist(turning_points[i]));
        total_cost += need_rot(dest_angle - angle);

        pos = turning_points[i];
        angle = dest_angle;
    }
    return total_cost;
}


// アイテムの取得順を求める
class OrderSolver
{
public:
    OrderSolver(const Stage& _stage, ActionSolver& _action_solver)
        : stage(_stage), items(stage.items()), n(items.count()), player(_stage.player()), action_solver(_action_solver)
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

    void set_order(int* o)
    {
        rep(i, n)
            order[i] = o[i];
    }

    void get_order(int* res) const
    {
        rep(i, n)
            res[i] = order[i];
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
    int calc_cost_simply() const
    {
        // 良いスコアが出るように調整した結果、この比率が良かった
        return 2 * calc_move_cost() + calc_rot_cost();
    }

    // まだ取得していないアイテムの中から最も距離が近いアイテムを順に取得していったときの順番を求める
    void greedy_near()
    {
        bool used[128] = {};

        Vec2 p = stage.player().pos();
        rep(oi, n)
        {
            int k = -1;
            float min_cost = 1e7;
            rep(i, n)
            {
                if (!used[i])
                {
                    float c = p.dist(items[i].pos());

                    // スタート地点からの移動のときにのみ、角度を考慮する
                    if (oi == 0)
                        c += abs(normalize_angle(get_angle(p, items[i].pos()))) * 9;

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
            used[k] = true;
        }
    }

    // 2-opt法 (2-optはTSPを解く時に用いられるヒューリスティックス)
    // スコア計算にはcalc_cost_simply()を使っている
    // calc_cost_simply()は正確性に欠けるが高速
    void improve_simply()
    {
        if (n < 3)
            return;

        int cur = calc_cost_simply();
        int ori[128];
        rep(i, n)
            ori[i] = order[i];
        for (;;)
        {
            bool updated = false;
            rep(b, n) rep(a, b - 1)
            {
                swap(order[a + 1], order[b]);
                reverse(order + a + 2, order + b);

                int next = calc_cost_simply();
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

    // 条件付きの2-opt
    // 厳密に2-optすると時間が足りないので、経路が交差しているものだけを評価する
    // スコア計算にはcalc_cost()を使っている
    // calc_cost()は実際にかかるターン数を正確に表すが低速
    void improve()
    {
        if (n < 3)
            return;

        int cur = calc_cost();
        int ori[128];
        rep(i, n)
            ori[i] = order[i];
        for (;;)
        {
            bool updated = false;
            for (int a = 0; a < n - 1; ++a)
            {
                for (int b = a + 2; b < n - 1; ++b)
                {
                    if (!intersect(pos[order[a]], pos[order[a + 1]], pos[order[b]], pos[order[b + 1]]))
                        continue;

                    swap(order[a + 1], order[b]);
                    reverse(order + a + 2, order + b);

                    int next = calc_cost();
                    if (next < cur)
                    {
                        updated = true;
                        cur = next;
                        rep(i, n)
                            ori[i] = order[i];
                        break;
                    }
                    else
                    {
                        rep(i, n)
                            order[i] = ori[i];
                    }
                }
            }
            if (!updated)
                break;
        }
    }

    int calc_cost()
    {
        Vec2 p[1024];
        int np = action_solver.solve(order, p);
        return solver::calc_cost(stage, p, np);
    }

private:
    const Stage& stage;
    const ItemCollection& items;
    const int n;
    const solver::Player player;
    ActionSolver& action_solver;

    int order[128];

    Vec2 pos[128];
    int mcost[128][128];
    float ang[128][128];
};


Answer answer;
int ans_i;
void solve(const Stage& stage)
{
    ActionSolver solver(stage);

    // アイテムの取得順を求める
    OrderSolver order_solver(stage, solver);
    int order[128];
    order_solver.greedy_near(); // 初期解
    order_solver.improve_simply(); // 高速なアルゴリズムで大まかに
    order_solver.improve(); // 低速だが正確なアルゴリズムで細かいスコア改善
    order_solver.get_order(order);

    // 回転する座標を求める
    Vec2 pos[1024];
    const int n = solver.solve(order, pos);

    solver::make_actions(stage, pos, n, answer);

    ans_i = 0;
}

hpc::Action get_next_action()
{
    const Action& t = answer.action(ans_i++);
    return hpc::Action(t.type(), t.value());
}

#ifdef LOCAL
void print_stage_info(const Stage& stage)
{
    const int move = answer.num_move();
    const int rot = answer.num_rot();
    fprintf(stderr, "stage score: %6.2f\n", answer.score(stage));
    fprintf(stderr, "move, rot: %5d, %5d\n", move, rot);
    debug(stage.items().count());
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
        debug(solver::stage_no);

        solver::solve(aStage);

#ifdef LOCAL
        solver::print_stage_info(aStage);
        cerr << endl;
#endif
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

        return solver::get_next_action();
    }
}

