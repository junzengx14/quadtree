import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

class Node2D:
    def __init__(self, box, depth, val):
        '''
        A 2D node in the two-dimensional quadtree
        '''
        self.box = box # leftbound, bottombound, rightbound, topbound
        self.depth = depth
        self.val = val
        self.children = [] # bottom_left, bottom_right, top_left, top_right
        self.check()

    def check(self):
        '''
        Check for some additional assumptions (TODO: Remove it)
        '''
        if (self.box[2] <= self.box[0] | self.box[3] <= self.box[1]):
            raise Exception("Invalid node, Current implementation assume that right/top bound should be larger than left/bottom bound")

    def split(self):
        '''
        Split the node with four children
        (l,t)----------(cx,t)----------(r,t)
                tl              tr
        (l,cy)---------(cx,cy)---------(r,cy)
                bl              br
        (l,b)----------(cx,b)----------(r,b)
        '''
        l, b, r, t = self.box
        # cx, cy are the coordinates of cell's center point 
        cx = l + (r - l) / 2
        cy = t + (b - t) / 2
        # generate four nodes connected to the current node
        tl = Node2D((l, cy, cx, t), self.depth + 1, self.val)
        tr = Node2D((cx, cy, r, t), self.depth + 1, self.val)
        bl = Node2D((l, b, cx, cy), self.depth + 1, self.val)
        br = Node2D((cx, b, r, cy), self.depth + 1, self.val)
        self.children = (bl, br, tl, tr)
    
    def get_children(self):
        '''
        Return children
        '''
        return self.children
    
    def get_child_bottomleft(self):
        '''
        Return the child at bottom left
        '''
        return self.children[0]
    
    def get_child_bottomright(self):
        '''
        Return the child at bottom right
        '''
        return self.children[1]

    def get_child_topleft(self):
        '''
        Return the child at top left
        '''
        return self.children[2]
    
    def get_child_topright(self):
        '''
        Return the child at top right
        '''
        return self.children[3]
    
    def get_bound_left(self):
        '''
        Return the left bound
        '''
        return self.box[0]
    
    def get_bound_bottom(self):
        '''
        Return the bottom bound
        '''
        return self.box[1]

    def get_bound_right(self):
        '''
        Return the right bound
        '''
        return self.box[2]
    
    def get_bound_top(self):
        '''
        Return the top bound
        '''
        return self.box[3]

    def get_child_from_point(self, point):
        '''
        Return the child which contains the given point
        '''
        if ( (0.5 * (self.get_bound_bottom() + self.get_bound_top()) < point[1]) & (point[1] < self.get_bound_top()) ):
            if ( (self.get_bound_left() < point[0]) & (point[0] < 0.5 * (self.get_bound_left() + self.get_bound_right())) ):
                return self.get_child_topleft()
            else:
                return self.get_child_topright()
        elif ( (self.get_bound_bottom() < point[1]) & (point[1] < 0.5 * (self.get_bound_bottom() + self.get_bound_top())) ):
            if ( (self.get_bound_left() < point[0]) & (point[0] < 0.5 * (self.get_bound_left() + self.get_bound_right())) ):
                return self.get_child_bottomleft()
            else:
                return self.get_child_bottomright()
        else:
            raise Exception("Cannot get child that contains the given point")

    def get_leaf_from_point(self, point):
        '''
        Return the leaf which contains the given point
        '''
        node = self
        while (node.is_leaf() != True):
            node = node.get_child_from_point(point)
        return node

    def get_leaves(self):
        '''
        Return the leaves
        ''' 
        leaves_list, node_list = [], []
        node_list.append(self)
        while (len(node_list) != 0):
            # BFS to get all leaves
            node = node_list.pop(0)
            if (node.is_leaf()):
                leaves_list.append(node)
            else:
                node_list.extend(node.get_children())
        return leaves_list

    def get_color(self):
        '''
        Return display from the node's value
        '''
        if (self.val == -1): # undefined
            return 'grey'
        elif (self.val == 0): # free
            return 'white'
        elif (self.val == 1): # occuped/obstacle
            return 'red'
        else:
            raise Exception("Undefined color for node")

    def get_rectangle(self):
        '''
        Return rectangle representation
        '''
        return Rectangle((self.box[0], self.box[1]), self.box[2]-self.box[0], self.box[3]-self.box[1], \
            edgecolor='black', facecolor = self.get_color())

    def generate_perfect_subquadtree(self, max_depth):
        '''
        Generate a perfect (sub)quadtree into a given depth from the node without children
        '''
        node_list = []
        node_list.append(self)
        while (len(node_list) != 0):
            # BFS for tree structure generation
            node = node_list.pop(0)
            if (node.depth < max_depth):
                node.split()
                leaves = node.get_children()
                node_list.extend(leaves)    

    def is_leaf(self):
        '''
        Return true or false if node is a leaf or not
        '''
        return (self.children == [])

    def is_contain_point(self, point):
        '''
        Return true or false if a point is in the node's representation (rectangle)
        '''
        if ((self.get_bound_left < point[0]) & (point[0] < self.get_bound_right)):
            if ((self.get_bound_bottom < point[1]) & (point[1] < self.get_bound_top)):
                return True
        return False

    def update_val(self, val):
        '''
        Update the node's value with a given value
        '''
        self.val = val # TODO: make it as continunous transition

class QuadTree2D():
    def __init__(self, root, max_depth):
        self.root = root
        self.max_depth = max_depth

    def explore_update(self, point, val):
        '''
        Given a point within the quadtree, explore/expand the subquadtrees and update the leaf value
        '''
        node = self.get_leaf_from_point(point)
        if (node.val == val):
            # same value, so that we don't need to expand
            pass
        else:
            while (node.depth < self.max_depth):
                # split the node
                node.split()
                # get the leaf which contains the given point
                node = node.get_child_from_point(point)
            # update the leaf's color
            node.update_val(val)

    def get_leaf_from_point(self, point):
        '''
        Return a leaf which contains the given point
        '''
        node = self.root
        while (node.children != []):
            node = node.get_child_from_point(point)
        return node

    def get_leaf_from_gridmap(self, index, depth):
        '''
        Return the leaf with the given x,y index from 2D grip map.

        Assumption:
            It assumes that the quadtree is fully expanded
        '''
        node = self.root
        x, y = int(index[0]), int(index[1])
        index_x_list, index_y_list = [], []
        for i in range(depth):
            index_x_list.append(x % 2)
            index_y_list.append(y % 2)
            x //= 2
            y //= 2
        for i in range(depth):
            # this trick comes from the children order
            node = node.get_children()[index_x_list[depth-i-1] + 2*index_y_list[depth-i-1]]
        return node

    def get_rectangle_from_root(self):
        '''
        Return the rectangle patch generated from the current node
        '''
        return self.root.get_rectangle()
        
    def get_rectangles_from_map(self):
        '''
        Return the list of rectangle patches generated from the current map
        '''
        rectangle_list = []
        node_list = self.root.get_leaves()
        for node in node_list:
            rectangle_list.append(node.get_rectangle())
        return rectangle_list

    def generate_perfect_quadtree(self):
        '''
        Generate the perfect quadtree from the root
        '''
        self.root.generate_perfect_subquadtree(self.max_depth)

    def plot(self):
        '''
        Visualize the quadtree
        '''
        fig = plt.figure()
        current_axes = fig.gca()
        rectangle_list = self.get_rectangles_from_map()
        for rec in rectangle_list:
            current_axes.add_patch(rec)
        current_axes.set_xlim(self.root.box[0], self.root.box[2])
        current_axes.set_ylim(self.root.box[1], self.root.box[3])
        plt.show()