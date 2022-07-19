# GPDynamics.py
# Marcio Gameiro
# 2022-07-19
# MIT LICENSE

import CMGDB
import DSGRN

import numpy as np

def Coordinates(box, phase_subdiv, lower_bounds, upper_bounds):
    # Returns coordinates of grid elements
    coord = [round((box[0] - lower_bounds[0]) * (2 ** phase_subdiv)/(upper_bounds[0] - lower_bounds[0]))]
    return coord

def get_intervals(morse_set):
    intervals = []
    x_min, x_max = 0, 0
    for box in sorted(morse_set):
        if x_min == x_max:
            # Initialize
            x_min, x_max = box[0:2]
            continue
        # If box glued to previous interval
        if x_max == box[0]:
            x_max = box[1]
        else:
            # Start new interval
            intervals.append([x_min, x_max])
            x_min, x_max = box[0:2]
    # Append final interval
    intervals.append([x_min, x_max])
    return intervals

def box_cover(domain_boxes, rect):
    """Find set of boxes intersecting the interval rect"""
    # Get indices of min and max vertices
    # Append -1 and n_verts to catch out of range cases
    n_verts = len(domain_boxes)
    v_min = max([-1] + [k for k, x in enumerate(domain_boxes) if x < rect[0]])
    v_max = min([n_verts] + [k for k, x in enumerate(domain_boxes) if x > rect[1]])
    # Return all indices between v_min and v_max
    # Return -1 and/or n_verts if rect is out of range
    return list(range(v_min, v_max))

def ComputeDomainGraph(phase_subdiv, lower_bounds, upper_bounds, g, confidence_level, L):
    # List of phase space boxes
    domain_boxes = np.linspace(lower_bounds[0], upper_bounds[0], 2**phase_subdiv + 1)
    # Center points of phase space boxes
    center_points = [sum(domain_boxes[k:k+2]) / 2 for k in range(len(domain_boxes) - 1)]
    # Map g evaluated at the center points
    g_center = np.array([(g([x])[0]) for x in center_points])
    # Standard deviation at the center points
    sigma_center = np.array([(g([x])[1]) for x in center_points])
    # Grid size h
    h = domain_boxes[1] - domain_boxes[0]
    # Get mean plus or minus confidence_level times std at center points
    g_l_center = g_center - confidence_level * sigma_center
    g_u_center = g_center + confidence_level * sigma_center
    # Number of points (edges vertices)
    num_points = len(domain_boxes)
    # Number of vertices (edges)
    num_verts = len(domain_boxes) - 1
    # Construct digraph
    domain_graph = DSGRN.Digraph()
    # Set number of vertices
    domain_graph.resize(num_verts)

    # Get MVM edges for odd cells
    for k1 in range(1, num_verts, 2):
        # Get confidence interval
        rect = [g_l_center[k1], g_u_center[k1]]
        # Get box cover of confidence interval
        b_cover = box_cover(domain_boxes, rect)
        # Discard out of range cases (-1 and num_points - 1)
        # List of grid boxes indices covering interval
        grid_cover = [v for v in b_cover if v not in [-1, num_points - 1]]
        if not grid_cover:
            print('Bad grid cover!')
        # Add grid cover to list of edges
        for index in grid_cover:
            domain_graph.add_edge(k1, index)

    # Get MVM edges for even cells
    for k1 in range(1, num_verts - 1, 2):
        # Lower point
        s1 = h - (g_l_center[k1 + 2] - g_l_center[k1]) / (2.0 * L)
        x = center_points[k1] + s1
        y = g_l_center[k1] - s1 * L
        # Make sure the rays intersect, that is, x is in
        # between the two consecutive odd center points
        if not (center_points[k1] <= x <= center_points[k1 + 2]):
            print('Rays do not intersect')
        y1 = y
        # Upper point
        s1 = h + (g_u_center[k1 + 2] - g_u_center[k1]) / (2.0 * L)
        x = center_points[k1] + s1
        y = g_u_center[k1] + s1 * L
        # Make sure the rays intersect, that is, x is in
        # between the two consecutive odd center points
        if not (center_points[k1] <= x <= center_points[k1 + 2]):
            print('Rays do not intersect')
        y2 = y
        # Get cone interval
        rect = [y1, y2]
        # Get box cover of confidence interval
        b_cover = box_cover(domain_boxes, rect)
        # Discard out of range cases (-1 and num_points - 1)
        # List of grid boxes indices covering interval
        grid_cover = [v for v in b_cover if v not in [-1, num_points - 1]]
        if not grid_cover:
            print('Bad grid cover!')
        # Add grid cover to list of edges
        for index in grid_cover:
            domain_graph.add_edge(k1 + 1, index)

    # Get MVM edges for first cell
    k1 = 1
    s = 3 * h / 2
    y1 = g_l_center[k1] - s * L
    y2 = g_u_center[k1] + s * L
    # Get cone interval
    rect = [y1, y2]
    # Get box cover of confidence interval
    b_cover = box_cover(domain_boxes, rect)
    # Discard out of range cases (-1 and num_points - 1)
    # List of grid boxes indices covering interval
    grid_cover = [v for v in b_cover if v not in [-1, num_points - 1]]
    if not grid_cover:
        print('Bad grid cover!')
    # Add grid cover to list of edges
    for index in grid_cover:
        domain_graph.add_edge(k1 - 1, index)

    return domain_graph, domain_boxes

def ConleyMorseGraph(phase_subdiv, lower_bounds, upper_bounds, g, confidence_level, L=None):
    """Compute Conley Morse graph for the GP map g"""

    if L == None:
        # Define box map
        def F(box):
            # End points of interval
            x1, x2 = box
            # Evaluate the map g at the end points
            y1, y2 = g([x1])[0], g([x2])[0]
            # Standard deviation at the end points
            sigma1, sigma2 = g([x1])[1], g([x2])[1]
            # Get end points of covering rectangle
            y_rect1 = min(y1 - confidence_level * sigma1, y2 - confidence_level * sigma2)
            y_rect2 = max(y1 + confidence_level * sigma1, y2 + confidence_level * sigma2)
            y_rect = [y_rect1, y_rect2]
            return y_rect
    else:
        # Compute domain graph using Lispschitz constant
        domain_graph, domain_boxes = ComputeDomainGraph(phase_subdiv, lower_bounds, upper_bounds, g, confidence_level, L)

        # Define a new box map from the multi-valued map constructed above.
        # The image of a box is just a rectangle containing all the boxes
        # in the image of the multi-valued map.

        # Grid size h
        h = domain_boxes[1] - domain_boxes[0]

        # Define box map
        def F(box):
            # Index of box (vertex u)
            u = Coordinates(box, phase_subdiv, lower_bounds, upper_bounds)[0]
            # Adjacencies of u
            adj_u = domain_graph.adjacencies(u)
            # Min and max vertices
            v_min, v_max = min(adj_u), max(adj_u)
            y_rect = [domain_boxes[v_min] + h / 2, domain_boxes[v_max + 1] - h / 2]
            return y_rect

    model = CMGDB.Model(phase_subdiv, lower_bounds, upper_bounds, F)
    morse_graph, map_graph = CMGDB.ComputeConleyMorseGraph(model)
    return morse_graph, map_graph
