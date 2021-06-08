import random, math
import numpy as np


class Constraint:
    def __init__(self, x, y, direction, b):
        self.x = x
        self.y = y
        self.direction = direction
        self.b = b


class Objective:
    def __init__(self, x, y, type):
        self.x = x
        self.y = y
        self.type = type
        self.b = 0


def intersection(a, b):
    global obj
    w = a.x * b.y - b.x * a.y
    wx = a.b * b.y - b.b * a.y
    wy = a.x * b.b - b.x * a.b
    if w != 0:
        return [wx / w, wy / w]
    else:
        return None


def vertex(points):
    maximum = [0, 0]
    global obj
    if obj.type == 'max':
        maximum[0] = points[0][0]
        maximum[1] = points[0][1]
        for p in points[1:]:
            if p[0] * obj.x + p[1] * obj.y > maximum[0] + maximum[1]:
                maximum[0] = p[0]
                maximum[1] = p[1]
        return maximum
    else:
        maximum[0] = points[0][0]
        maximum[1] = points[0][1]
        for p in points[1:]:
            if p[0] * obj.x + p[1] * obj.y < maximum[0] + maximum[1]:
                maximum[0] = p[0]
                maximum[1] = p[1]
        return maximum


def check_if_in_halfplane(c, v):
    if c.direction == 'g':
        if v[0] * c.x + v[1] * c.y >= c.b:
            return True
        else:
            return False
    else:
        if v[0] * c.x + v[1] * c.y <= c.b:
            return True
        else:
            return False


def pre_process():
    global obj, constraints
    v1 = [obj.x, obj.y]
    v2 = [constraints[0].x, constraints[0].y]
    dp = np.dot(v1 / np.linalg.norm(v1), v2 / np.linalg.norm(v2))
    minimum = np.arccos(dp)
    if minimum > math.pi:
        minimum -= math.pi
    s = constraints[0]
    dps = dp
    for c in constraints[1:]:  # 1.
        v1 = [obj.x, obj.y]
        v2 = [c.x, c.y]
        dp = np.dot(v1 / np.linalg.norm(v1), v2 / np.linalg.norm(v2))
        min2 = np.arccos(dp)
        if min2 > math.pi:
            min2 -= math.pi
        if minimum > min2:
            minimum = np.arccos(dp)
            s = c
            dps = dp
    if s != constraints[0]:
        i = constraints.index(s)
        constraints[0], constraints[i] = constraints[i], constraints[0]
    for c in constraints[1:]:  # 2.
        v1 = [obj.x, obj.y]
        v2 = [c.x, c.y]
        dp = np.dot(v1 / np.linalg.norm(v1), v2 / np.linalg.norm(v2))
        if abs(dps - dp) < 0.0001 and not check_if_in_halfplane(s, v2):
            return "infeasible"
    for c in constraints[1:]:  # 3,4
        p = intersection(s, c)
        if obj.type == 'max':
            p[0] += obj.x
            p[1] += obj.y
            if check_if_in_halfplane(c, p):
                continue
            else:
                i = constraints.index(c)
                constraints[1], constraints[i] = constraints[i], constraints[1]
                return "ok"
        else:
            p[0] -= obj.x
            p[1] -= obj.y
            if check_if_in_halfplane(c, p):
                continue
            else:
                i = constraints.index(c)
                constraints[1], constraints[i] = constraints[i], constraints[1]
                return "ok"
    return "unbounded"


def seidel(constraints):
    global obj
    points = []
    z = pre_process()
    if z == "ok":
        i = 2
        v = intersection(constraints[0], constraints[1])
        tmp = constraints[2:]
        random.shuffle(tmp)
        constraints[2:] = tmp
        for c in constraints[2:]:
            i += 1
            if not check_if_in_halfplane(c, v):
                for l in constraints[:i]:
                    if l != c and intersection(c, l) is not None:
                        points.append(intersection(c, l))
                points2 = points[:]
                for p in points2:
                    for h in constraints[:i]:
                        if not check_if_in_halfplane(h, p):
                            if p in points:
                                points.remove(p)
                if not points:
                    return "infeasible"
                v = vertex(points)
                points.clear()
        return v
    else:
        return z


constraints = []
obj = Objective(3, 4, 'max')

constraints.append(Constraint(2, 3, 'l', 4))
constraints.append(Constraint(2, 5, 'l', 5))
constraints.append(Constraint(2, 1, 'l', 3))
constraints.append(Constraint(3, 4, 'l', 2))
constraints.append(Constraint(-1, 0, 'l', 0))
constraints.append(Constraint(0, -1, 'l', 0))


print(seidel(constraints))
