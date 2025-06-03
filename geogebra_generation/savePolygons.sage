from numpy import random
from sage.geometry.polyhedron.plot import cyclic_sort_vertices_2d;

def construct_polygons(list_of_vertices, bound_box = Polyhedron([(10000,10000),(10000,-10000),(-10000,-10000),(-10000,10000)])):
    colors = ["blue", "green", "orange", "purple", "cyan"]
    poly = polygon(bound_box, alpha=0)

    for points in list_of_vertices:
        if ("inf", "inf") in points:
            continue

        points = [(Rational(float(p[0])), Rational(float(p[1]))) for p in points]
        P = Polyhedron(points).intersection(bound_box)
        vs = cyclic_sort_vertices_2d(P.vertices())
        if len(vs) > 2:
            poly = poly + polygon(vs, color=list(random.rand(3,)), dpi=4000)
    
    return poly


if __name__ == "__main__":
    all_points = []
    with open("save_regions_test.txt", "r") as p:
        for i, p in enumerate(p.readlines()):
            data = p.strip("\n").split(":")
            if len(data[1]) == 0:
                continue
            points = [tuple(p.strip("()").replace(" ", "").split(",")) for p in data[1].split("), (")]
            all_points.append(points)


    poly = construct_polygons(all_points)
            # poly += point((56/5,-67/10), color="black",size=1,zorder=100)

            # for p in [(-50,103/10),(-50,323/30),(-50,173/15),(-50,361/30),(-50,182/15),(-50,188/15),(-50,193/15),(-50,403/30),(-50,419/30),(-50,147/10),(-50,16),(-50,81/5),(-50,88/5),(-50,187/10),(50,-173/26),(50,-7),(50,-73/10),(50,-39/5),(50,-42/5),(50,-91/10),(50,-93/10),(50,-47/5),(50,-101/10),(50,-56/5),(50,-121/10),(50,-253/10),(50,-128/5),(50,-521/20),(50,-261/10),(50,-529/20),(50,-53/2),(50,-136/5),(50,-551/20),(50,-761/25),(50,-169/5),(50,-819/20),(-50,171/10)]:
            #     poly += point(p, color="black",size=1,zorder=100)
    poly.save("polygraph.png")