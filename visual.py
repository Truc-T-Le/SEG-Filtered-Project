import math
import random
import os
import argparse
from PIL import Image, ImageDraw, ImageFont
import assets.circlify as circ
import pandas as pd
import csv

# GLOBAL TO TRACK FILES
fid_to_location = {}

"""
# UTIL FUNCTIONS
"""
def parse_args():
    parser = argparse.ArgumentParser(description="Analysis program.")
    parser.add_argument("input", help="Directory containing all input files")
    parser.add_argument("--out", help="Output image name", default="output.png")
    args = parser.parse_args()

    return args

"""
# GROUP FUNCTIONS
"""
def get_group_name(pair):
    if cation_pi_check(pair):
        return "Cation-Pi"

    groups = []

    for letter in pair:
        if hydro_check(letter):
            groups.append("Hydro")
        elif aromatic_check(letter):
            groups.append("Aromatic")
        elif polar_check(letter):
            groups.append("Polar")
        elif charge_check(letter):
            groups.append("Charge")
        elif unique_check(letter):
            groups.append("Unique")

    return "-".join(groups)

def cation_pi_check(pair):
    # Special check
    letters = "KRYF"
    for pair_letter in pair:
        if pair_letter in letters:
            return True
    return False

def hydro_check(letter):
    letters = "AILMV"
    return letter in letters

def aromatic_check(letter):
    letters = "FYW"
    return letter in letters

def polar_check(letter):
    letters = "QNSTC"
    return letter in letters

def charge_check(letter):
    letters = "RKEDH"
    return letter in letters

def unique_check(letter):
    letters = "GP"
    return letter in letters

"""
# UTIL FUNCTIONS
"""

def circle_pos(pos, radius, angle):
    x = pos[0] + radius * math.cos(angle)
    y = pos[1] + radius * math.sin(angle)
    return (x, y)

def rand_angle():
    return random.random() * (2 * math.pi)

class Circle:
    def __init__(self, position, radius):
        self.position = position
        self.radius = radius

    def single_circle_overlapping(self, circle):
        dist = math.dist(self.position, circle.position)
        return dist < self.radius + circle.radius

    def circles_not_overlapping(self, circles):
        for circle in circles:
            if self.single_circle_overlapping(circle):
                return False
        return True

    def circles_not_in_canvas(self, border):
        if self.position[0] - self.radius < border[0]:
            return True
        if self.position[0] + self.radius > border[1]:
            return True
        if self.position[1] - self.radius < border[2]:
            return True
        if self.position[1] + self.radius > border[3]:
            return True
        return False

def repeat_draw_until_success(randomize_fn, draw_fn, circles_drawn, border_check_args=None):
    for iteration in range(100):
        circle = randomize_fn()
        if circle.circles_not_overlapping(circles_drawn):
            if border_check_args != None:
                if not circle.circles_not_in_canvas(border_check_args):
                    break
            else:
                break


    draw_fn(circle)
    return circles_drawn + [circle] 

"""
# DATA STRUCT
"""
class SingleNode:
    def __init__(self, indexed_id):
        self.fid = indexed_id.split(".")[0]
        self.index = indexed_id.split(".")[1]
        self.percentage = None
        self.normalized_percentage = None

    def calc_percentage(self, pair, culled_str):
        # Make sure they're both lower case
        culled_str = culled_str.lower()
        pair = pair.lower()

        occurences = culled_str.count(pair[0]) + culled_str.count(pair[1])
        self.percentage = round((occurences / len(culled_str)) * 100)

    def get_radius(self):
        return 5

    def normalize_percentage(self, min_percentage, max_percentage):
        range_percentage = max_percentage - min_percentage
        if range_percentage != 0:
            self.normalized_percentage = (self.percentage - min_percentage) / range_percentage
        else:
            self.normalized_percentage = 0

    def draw(self, draw_ctx, center, network_type):
        # Log position to draw lines
        if self.fid not in fid_to_location:
            fid_to_location[self.fid] = []
        fid_to_location[self.fid].append(center)

        single_node_radius = self.get_radius()
        position_a = (center[0] - single_node_radius, center[1] - single_node_radius)
        position_b = (center[0] + single_node_radius, center[1] + single_node_radius)
        draw_ctx.ellipse([position_a, position_b], fill=(40, 82, 122, 80))

class ClusterNode:
    def __init__(self, pair):
        self.pair = pair

        # Map from indexed_id ("15.1") to the Single Nodes
        self.single_nodes = {}
    def insert(self, indexed_id, culled_str):
        # Indexed id must be unique
        assert(indexed_id not in self.single_nodes)

        self.single_nodes[indexed_id] = SingleNode(indexed_id)
        self.single_nodes[indexed_id].calc_percentage(self.pair, culled_str)

    def get_radius(self):
        # NOTE: Log used so that radius doesn't grow too large
        return 30 * math.log(1 + len(self.single_nodes))

    def draw(self, draw_ctx, center, network_type):
        # Sort single_nodes by percentage
        sorted_items = sorted(self.single_nodes.items(), key=lambda x: x[1].percentage)
        for node in self.single_nodes.values():
            node.normalize_percentage(sorted_items[0][1].percentage, sorted_items[-1][1].percentage)
        
        # Draw cluster circle
        cluster_radius = self.get_radius()
        position_a = (center[0] - cluster_radius, center[1] - cluster_radius)
        position_b = (center[0] + cluster_radius, center[1] + cluster_radius)
        if network_type == "gcs":
            font = ImageFont.truetype("assets/Roboto-Light.ttf", 25)
            text_pos = (center[0] - 20, center[1] + cluster_radius)
            draw_ctx.text(text_pos, self.pair, font=font, fill="BLACK")
            draw_ctx.ellipse([position_a, position_b], fill=(138, 196, 208, 123))
        elif network_type == "gs":
            pass



        # Draw nodes within cluster
        circles_drawn = []
        for indexed_id, single_node in sorted_items:
            def randomize_fn():
                single_node_radius = single_node.get_radius()
                radius = (cluster_radius - single_node_radius) * single_node.normalized_percentage 
                angle = rand_angle()
                cur_position = circle_pos(center, radius, angle)

                return Circle(cur_position, single_node_radius)

            def draw_fn(circle):
                single_node.draw(draw_ctx, circle.position, network_type)

            circles_drawn = repeat_draw_until_success(randomize_fn, draw_fn, circles_drawn)

class GroupNode:
    def __init__(self, group_name):
        self.group_name = group_name

        # Map from cluster name ("AB") to Cluster Nodes
        self.cluster_nodes = {}

    def insert(self, pair, indexed_id, culled_str):
        if pair not in self.cluster_nodes:
            self.cluster_nodes[pair] = ClusterNode(pair)

        self.cluster_nodes[pair].insert(indexed_id, culled_str)

    def get_radius(self):
        return sum(list(map(lambda x: x.get_radius() * 1.3, self.cluster_nodes.values())))

    def draw(self, draw_ctx, center, network_type):
        # Draw group circle
        group_radius = self.get_radius()
        position_a = (center[0] - group_radius, center[1] - group_radius)
        position_b = (center[0] + group_radius, center[1] + group_radius)
        draw_ctx.ellipse([position_a, position_b], fill=(244, 209, 96, 125))

        # Draw label
        font = ImageFont.truetype("assets/Roboto-Light.ttf", 30)
        text_pos = (center[0] - 80, center[1] - group_radius - 35)
        draw_ctx.text(text_pos, self.group_name, font=font, fill="BLACK")
        #print(self.group_name)

        # Draw nodes within group
        circles_drawn = []
        for key in self.cluster_nodes:
            def randomize_fn():
                cluster_radius = self.cluster_nodes[key].get_radius()
                radius = (group_radius - cluster_radius) * random.random()
                angle = rand_angle()
                cur_position = circle_pos(center, radius, angle)
                return Circle(cur_position, cluster_radius)

            def draw_fn(circle):
                self.cluster_nodes[key].draw(draw_ctx, circle.position, network_type)

            circles_drawn = repeat_draw_until_success(randomize_fn, draw_fn, circles_drawn)

class Network:
    def __init__(self):
        # Map from group name ("Hydro-Aro") to Group Nodes
        self.group_nodes = {}

    def insert(self, indexed_id, pair, culled_str):
        pair_sorted = "".join(sorted(pair))
        group_name = get_group_name(pair_sorted)

        if group_name not in self.group_nodes:
            self.group_nodes[group_name] = GroupNode(group_name)

        self.group_nodes[group_name].insert(pair_sorted, indexed_id, culled_str)


    def draw(self, network_type):
        # Sizing for the png outputs
        width = height = int(len(self.group_nodes) ** 1.30) * 120
        margin = width / 4

        width_start = margin / 2
        width_end = width - margin / 2
        width_span = width - margin
        height_start = margin / 2
        height_end = height - margin / 2
        height_span = height - margin
        borders = [width_start, width_end, height_start, height_end]

        image = Image.new("RGB", (width, height), "white")
        draw_ctx = ImageDraw.Draw(image, "RGBA")

        circles_drawn = []
        for key in self.group_nodes:
            def randomize_fn():
                radius = self.group_nodes[key].get_radius()
                center = (random.random() * width_span + width_start, random.random() * height_span + height_start)
                return Circle(center, radius)

            def draw_fn(circle):
                self.group_nodes[key].draw(draw_ctx, circle.position, network_type)

            circles_drawn = repeat_draw_until_success(randomize_fn, draw_fn, circles_drawn, borders)

        # Draw all lines
        global fid_to_location
        #print(fid_to_location)
        for dics, (names , centers) in enumerate(fid_to_location.items()):
            if len(centers) <= 1:
                continue
            for i in range(len(centers)):
                for j in range(i, len(centers)):
                    draw_ctx.line([centers[i], centers[j]], fill=(0, 0, 0, 60), width=3)
                # font = ImageFont.truetype("assets/Roboto-Light.ttf", 10)
                # text_pos = ([centers[i][0], centers[i][1]])
                # draw_ctx.text(text_pos, names, font=font, fill="BLACK")


        
        # Reset lines
        fid_to_location = {}

        return image

def net_counts(name_list, loc_list, net_dict):
    for idx, name in enumerate(name_list):
        if name in net_dict:
            net_dict[name].append(loc_list[idx]) # append the element
        else:
            net_dict[name] = [loc_list[idx]]
    # for keys, values in net_dict.items():
    #     if len(values) > 1:
    #         key = list(net_dict.keys())[list(net_dict.values()).index(values)]
    #         #print("seq")
    #         print(key, values)
    #     else:
    #         continue

    df = pd.DataFrame(net_dict.items(), columns = ["Sequence Name", "Classification"])
    df['Classification'] = df['Classification'].astype(str).str[1:-1]
    df.to_csv("Networks.csv", index=False)
    net_dict = {}

def net_counts_2(name_list, loc_list, net_dict):
    for idx, loc in enumerate(loc_list):
        if loc in net_dict:
            net_dict[loc].append(name_list[idx])
        else:
            net_dict[loc] = [name_list[idx]]
    df = pd.DataFrame(net_dict.items(), columns = ["Classification", "Sequence Name"])
    df['Sequence Name'] = df['Sequence Name'].astype(str).str[1:-1]
    df.to_csv("Networks_seq.csv", index=False)
    net_dict = {}


if __name__ == "__main__":
    args = parse_args()
    network = Network()
    cluster_input_dir = os.path.join(args.input, "clusters")
    ll = os.listdir(cluster_input_dir)
    clusters, groups = {}, {}
    net_dict = {}
    names, loc = [], []
    names2, loc2 = [], []
    for cluster_file in ll:
        cluster_name = cluster_file.split(".")[0]
        group_name = get_group_name(cluster_name)
        try:
            with open(os.path.join(cluster_input_dir, cluster_file), 'r') as f:
                cluster_lines = f.readlines()
        except:
            #print("Invalid file format found: {}".format(cluster_input_dir))
            continue
        groups[group_name] = 0
        clusters[cluster_name] = 0
        for cluster_line in cluster_lines[1:]:
            gc, gc2 = [], []
            groups[group_name] += 1
            clusters[cluster_name] += 1
            spl = cluster_line.split(",")
            network.insert(spl[0], cluster_name, spl[1])

            names.append(spl[0].split(".", 1)[0])
            gc.append(group_name)
            gc.append(cluster_name)
            loc.append(gc)
            names2.append(spl[0])
            #gc2.append(group_name + "__" + cluster_name)
            #gc2.append(cluster_name)
            loc2.append(group_name + " : " + cluster_name)



    net_counts(names, loc, net_dict)
    net_dict = {}
    net_counts_2(names2, loc2, net_dict)



    gcs_image = network.draw(network_type="gcs")
    gcs_image.save("gcs-"+ args.out, "PNG", dpi=(100000,100000))

    gc_image = network.draw(network_type="gs")
    gc_image.save("gs-"+ args.out, "PNG", dpi=(100000,100000))

    # circlify_data = []
    # for group_key, group_value in groups.items():
    #     circlify_data.append({'id': group_key, 'datum': group_value, 'children': []})
    #     for cluster_key, cluster_value in clusters.items():
    #         group_name = get_group_name(cluster_key)
    #         if group_name != group_key:
    #             continue
    #         circlify_data[-1]['children'].append({'id': cluster_key, 'datum': cluster_value})

    # data = [
    #     {'datum': 19, 'children': [0.5, 12]},
    #     {'datum': 5, 'children': [0.5]},
    #     {'datum': 5, 'children': [0.5]},
    # ]

    # circles = circ.circlify(circlify_data)
    # fname = "gc-"+ args.out
    # circ.bubbles(circles, fname)

