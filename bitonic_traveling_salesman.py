import math
from matplotlib import pyplot as plt

pts = [(1, 7), (2, 1), (3, 4), (6, 5), (7, 2), (8, 6), (9, 3)]

def euclidean_distance(point1, point2):
    return math.sqrt((point1[0] - point2[0]) ** 2 + (point1[1] - point2[1]) ** 2)

def bt4(points):
    bitonic_distance = [[0 for x in range(len(points))] for x in range(len(points))]
    tracker = [[[] for x in range(len(points))] for x in range(len(points))]
    lines = [[[] for x in range(len(points))] for x in range(len(points))]
    number_of_points = len(points)
    inner_loop_count = 0
    outer_loop_count = 0
    for i in range(1, number_of_points):
        inner_loop_count += 1
        bitonic_distance[0][i] = bitonic_distance[0][i - 1] + euclidean_distance(points[i - 1], points[i])
        tracker[0][i] = [0, i-1]
        lines[0][i] = [i-1, i]
    for i in range(1, number_of_points):
        for j in range(i, number_of_points):
            minimum = math.inf
            if (i == j) or (i == j - 1):
                for k in range(i):
                    outer_loop_count += 1
                    value_to_test = bitonic_distance[k][i] + euclidean_distance(points[k], points[j])
                    if value_to_test < minimum:
                        minimum = value_to_test
                        tracker[i][j] = [k, i]
                        lines[i][j] = [k, j]
                bitonic_distance[i][j] = minimum
            else:
                inner_loop_count += 1
                bitonic_distance[i][j] = bitonic_distance[i][j - 1] + euclidean_distance(points[j - 1], points[j])
                tracker[i][j] = [i, j-1]
                lines[i][j] = [j-1, j]
    print(inner_loop_count, outer_loop_count)
    return bitonic_distance, tracker, lines

def add_line_by_pt(pt1, pt2):
    plt.plot((pt1[0], pt2[0]), (pt1[1], pt2[1]), 'ro-')
    
def generate_plot(points, tracker, lines):
    n = len(points)
    start_end_point_index = lines[n-1][n-1]
    pt1 = start_end_point_index[0]
    pt2 = start_end_point_index[1]
    coord = tracker[n-1][n-1]
    add_line_by_pt(points[pt1], points[pt2])
    while n > 1:
        start_end_point_index = lines[coord[0]][coord[1]]
        # print(start_end_point_index)
        pt1 = start_end_point_index[0]
        pt2 = start_end_point_index[1]
        coord = tracker[coord[0]][coord[1]]
        add_line_by_pt(points[pt1], points[pt2])
        n -=1
    plt.show()
    

def main():
    distance, tracker, lines = bt4(pts)
    generate_plot(pts, tracker, lines)

main()