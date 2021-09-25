//
// Created by 159-mrv on 9/7/18.
//

#include <core/vertex/vertex.h>

ostream& operator<<(ostream& out, const ClearVertex& vertex) {
    out << vertex.index << " " << setprecision(30) << vertex.x << " " << vertex.y << " " << vertex.z;
    return out;
}