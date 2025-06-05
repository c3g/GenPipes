#!/usr/bin/env python3

import json
import os
import sys

import questionary 
from jinja2 import Environment 

def load_guide (file_path):
    """
    Load wizard JSON files by filename: {general, deployment, pipeline, protocol, command, step}_guide.json
    """
    full_file_path = os.path.join(os.path.dirname(__file__), "wizard_json", file_path)
    with open(full_file_path) as file:
        return json.load(file)
    
class Wizard:

    def __init__ (self, start_file):
        #store variables: pipeline_name, protocol_name, {r,j,c,o,d,p,s,g}_command, scheduler_server_name, 
        #server_in, path_custom_ini, final_command, step_range
        self.variables = {}

        #for rendering
        self.env = Environment()

        #load initial tree
        self.current_guide = load_guide(start_file)

        #for json file switching
        self.current_file = start_file

        #determines starting node based on node id
        self.current_node_id = self.current_guide["_meta"]["entry_point"]

    def apply_variables (self, message):
        """
        Fill the placeholder {{...}} in the message with the data from the current variable
        """
        return self.env.from_string(message).render(**self.variables)
    
    def find_node (self, node_id):
        """
        Find and return the node dictionary with the given id of the node
        """
        for n in self.current_guide["nodes"]:
            if n["id"] == node_id:
                return n
        raise RuntimeError(f"Node '{node_id}' not found in {self.current_file}")
    
    def goto(self, next_node):
        """
        Move to the next node in the tree 
        """
        # when next node is in the current JSON file
        if isinstance(next_node, str):
            self.current_node_id = next_node
        else:
            #when next node is in another JSON file
            self.current_file = next_node["external"]
            self.current_guide = load_guide(self.current_file)
            self.current_node_id = next_node["entryPoint"]

    def tree_traversal(self):
        """
        Traverse through the JSON files to prompt the user questions
        """

        #Keep asking questions/traversing the tree until reach the end
        while True:
            node = self.find_node(self.current_node_id)
            node_type = node.get("type")

            #Confirm: yes/no questions
            if node_type == "confirm":
                answer = questionary.confirm(self.apply_variables(node["question"])).ask()
                chosen = "Yes" if answer else "No"
                next_info = None
                for option in node["options"]:
                    if option["label"] == chosen:
                        next_info = option["next"]
                        break
                self.goto(next_info)

            #Selection: single-select question from list 
            elif node_type == "selection":
                #choices
                if "choices" in node:
                    labels = [c["label"] for c in node["choices"]]
                    choice = questionary.select(self.apply_variables(node["question"]), choices = labels).ask()
                    next_node = next(c for c in node["choices"] if c["label"] == choice)
                    self.goto(next_node["node"])

                #choices_cases
                else:
                    for case_block in node ["choices_cases"]:
                        variable_name, value = next(iter(case_block["when"]["equals"].items()))
                        if self.variables.get(variable_name) == value:
                            labels = [current_choice["label"] for current_choice in case_block["choices"]]
                            choice = questionary.select(self.apply_variables(node["question"]), choices = labels).ask()
                            next_node = next(current_choice for current_choice in case_block["choices"] if current_choice["label"] == choice)
                            self.goto(next_node["node"])

            #Set_variable: set and store variable 
            elif node_type == "set_variable":
                variable = node["variable"]
                raw_value = node["value"]
                updated_value = self.apply_variables(raw_value)
                self.variables[variable] = updated_value
                self.goto(node["next"])

            #Message: output message for user, if no next node then end of wizard
            elif node_type == "message":
                print (self.apply_variables(node["message"]))
                next_node = node.get("next")
                #end of wizard 
                if not next_node:
                    break
                self.goto(next_node)

            #Switch: determine the next node based on cases (e.g pipeline/protocol name)
            elif node_type == "switch":
                variable = node["variable"]
                value = self.variables.get(variable)
                cases = node["cases"]
                if value in cases:
                    self.goto(cases[value]["node"])
                else:
                    print(f"[ERROR] No matching case for {variable} ='{value}' at node {node['id']}")
                    sys.exit(1)

            #Input: Prompt the user for input and store it as a variable 
            elif node_type == "input":
                variable = node["variable"]
                input = questionary.text(self.apply_variables(node["prompt"])).ask()
                self.variables[variable] = input
                self.goto(node["next"])

            else:
                print(f"[ERROR] Unknown node type: {node_type} in {self.current_file}")
                sys.exit(1)

    #for testing
    def main():
        start_json_file = "general_guide.json"
        start = Wizard(start_json_file)
        start.run()

    if __name__ == "__main__":
        main()