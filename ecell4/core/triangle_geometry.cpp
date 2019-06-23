#include "triangle_geometry.hpp"
#include <cmath>

namespace ecell4
{

Real3 closest_point(const Real3& pos, const Triangle& tri)
{
    // this implementation is based on "Real-Time Collision Detection"
    // by Christer Ericson,
    // published by Morgan Kaufmann Publishers, (c) 2005 Elsevier Inc.
    // pp.141-142

    const Real3 a = tri.vertices()[0];
    const Real3 b = tri.vertices()[1];
    const Real3 c = tri.vertices()[2];

    const Real3 ab = b - a;
    const Real3 ac = c - a;
    const Real3 ap = pos - a;
    const Real  d1 = dot_product(ab, ap);
    const Real  d2 = dot_product(ac, ap);
    if (d1 <= 0.0 && d2 <= 0.0)
        return a;

    const Real3 bp = pos - b;
    const Real  d3 = dot_product(ab, bp);
    const Real  d4 = dot_product(ac, bp);
    if (d3 >= 0.0 && d4 <= d3)
        return b;

    const Real vc = d1*d4 - d3*d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0)
    {
        const Real v = d1 / (d1 - d3);
        return a + ab * v;
    }

    const Real3 cp = pos - c;
    const Real  d5 = dot_product(ab, cp);
    const Real  d6 = dot_product(ac, cp);
    if (d6 >= 0.0 && d5 <= d6)
        return c;

    const Real vb = d5*d2 - d1*d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0)
    {
        const Real w = d2 / (d2 - d6);
        return a + ac * w;
    }

    const Real va = d3*d6 - d5*d4;
    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0)
    {
        const Real w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return b + (c - b) * w;
    }

    const Real denom = 1.0 / (va + vb + vc);
    const Real v     = vb * denom;
    const Real w     = vc * denom;
    return a + ab * v + ac * w;
}

Real3 closest_point(const Real3& pos, const TriangleView& tri)
{
    // this implementation is based on "Real-Time Collision Detection"
    // by Christer Ericson,
    // published by Morgan Kaufmann Publishers, (c) 2005 Elsevier Inc.
    // pp.141-142

    const Real3 a = tri.vertices(0);
    const Real3 b = tri.vertices(1);
    const Real3 c = tri.vertices(2);

    const Real3 ab = b - a;
    const Real3 ac = c - a;
    const Real3 ap = pos - a;
    const Real  d1 = dot_product(ab, ap);
    const Real  d2 = dot_product(ac, ap);
    if (d1 <= 0.0 && d2 <= 0.0)
        return a;

    const Real3 bp = pos - b;
    const Real  d3 = dot_product(ab, bp);
    const Real  d4 = dot_product(ac, bp);
    if (d3 >= 0.0 && d4 <= d3)
        return b;

    const Real vc = d1*d4 - d3*d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0)
    {
        const Real v = d1 / (d1 - d3);
        return a + ab * v;
    }

    const Real3 cp = pos - c;
    const Real  d5 = dot_product(ab, cp);
    const Real  d6 = dot_product(ac, cp);
    if (d6 >= 0.0 && d5 <= d6)
        return c;

    const Real vb = d5*d2 - d1*d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0)
    {
        const Real w = d2 / (d2 - d6);
        return a + ac * w;
    }

    const Real va = d3*d6 - d5*d4;
    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0)
    {
        const Real w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return b + (c - b) * w;
    }

    const Real denom = 1.0 / (va + vb + vc);
    const Real v     = vb * denom;
    const Real w     = vc * denom;
    return a + ab * v + ac * w;
}

Real3 closest_point(const Real3& pos, const TriangleConstView& tri)
{
    // this implementation is based on "Real-Time Collision Detection"
    // by Christer Ericson,
    // published by Morgan Kaufmann Publishers, (c) 2005 Elsevier Inc.
    // pp.141-142

    const Real3 a = tri.vertices(0);
    const Real3 b = tri.vertices(1);
    const Real3 c = tri.vertices(2);

    const Real3 ab = b - a;
    const Real3 ac = c - a;
    const Real3 ap = pos - a;
    const Real  d1 = dot_product(ab, ap);
    const Real  d2 = dot_product(ac, ap);
    if (d1 <= 0.0 && d2 <= 0.0)
        return a;

    const Real3 bp = pos - b;
    const Real  d3 = dot_product(ab, bp);
    const Real  d4 = dot_product(ac, bp);
    if (d3 >= 0.0 && d4 <= d3)
        return b;

    const Real vc = d1*d4 - d3*d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0)
    {
        const Real v = d1 / (d1 - d3);
        return a + ab * v;
    }

    const Real3 cp = pos - c;
    const Real  d5 = dot_product(ab, cp);
    const Real  d6 = dot_product(ac, cp);
    if (d6 >= 0.0 && d5 <= d6)
        return c;

    const Real vb = d5*d2 - d1*d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0)
    {
        const Real w = d2 / (d2 - d6);
        return a + ac * w;
    }

    const Real va = d3*d6 - d5*d4;
    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0)
    {
        const Real w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return b + (c - b) * w;
    }

    const Real denom = 1.0 / (va + vb + vc);
    const Real v     = vb * denom;
    const Real w     = vc * denom;
    return a + ab * v + ac * w;
}

} // ecell4
