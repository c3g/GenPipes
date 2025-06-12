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
        self.variables = {}
        self.env = Environment()
        self.current_guide = load_guide(start_file)
        self.current_file = start_file
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

                if variable in ("r_command", "path_custom_ini", "g_command", "d_command", "p_command"):
                    self.fix_filenames()

                if "o_command" in variable:
                    self.fix_filenames()

                updated_value = self.apply_variables(raw_value)

                #Clean up extra spaces in final command
                if variable == "final_command":
                    updated_value = " ".join(updated_value.split())

                #updated_value = self.apply_variables(raw_value)
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
                while True:
                    input = questionary.text(self.apply_variables(node["prompt"])).ask()

                    self.variables[variable] = input

                    if variable == "step_range":
                        self.variables[variable] = input
                        if not self.valid_step_range():
                            continue
                    break
                
                self.goto(node["next"])

            else:
                print(f"[ERROR] Unknown node type: {node_type} in {self.current_file}")
                sys.exit(1)
    
    def fix_filenames(self):
        """
        Handle cases where user includes/doesn't include .txt/ini/sh to their input
        """
        readset_filename = self.variables.get("raw_readset_filename", "").strip()
        if readset_filename and not readset_filename.endswith(".txt"):
            readset_filename += ".txt"
        self.variables["raw_readset_filename"] = readset_filename

        design_filename = self.variables.get("design_file_name", "").strip()
        if design_filename and not design_filename.endswith(".txt"):
            design_filename += ".txt"
        self.variables["design_file_name"] = design_filename

        pair_filename = self.variables.get("pair_file_name", "").strip()
        if pair_filename and not pair_filename.endswith(".txt"):
            pair_filename += ".txt"
        self.variables["pair_file_name"] = pair_filename


        path_custom_ini = self.variables.get("raw_path_custom_ini", "").strip()
        if path_custom_ini and not path_custom_ini.endswith(".ini"):
            path_custom_ini += ".ini"
        self.variables["raw_path_custom_ini"] = path_custom_ini

        genpipes_filename = self.variables.get("g_filename", "").strip()
        if genpipes_filename and not genpipes_filename.endswith(".sh"):
            genpipes_filename += ".sh"
        self.variables["g_filename"] = genpipes_filename

        #If pipeline has no protocol--> dont want -t in the final command
        t_command = self.variables.get("t_command", "").strip()
        if not self.variables.get("protocol_name"):
            t_command = ""
        self.variables["t_command"] = t_command

        #If user skips o command --> dont want -o in the final command
        o_command = self.variables.get("o_command", "").strip()
        if not self.variables.get("directory_name"):
            o_command = ""
        self.variables["o_command"] = o_command

    def valid_step_range(self):
        """
        Ensure that inputted step range is valid
        """
        pipeline = self.variables.get("pipeline_name")
        protocol = self.variables.get("protocol_name", "no_protocol")
        step_range = self.variables.get("step_range", "")

        valid_steps = {
            "ampliconseq": {"no_protocol": (1,8)},
            "chipseq": {
                "chipseq": (1,23),
                "atacseq": (1,24)
            },
            "covseq": {"no_protocol": (1, 21)},
            "dnaseq": {
                "germline_snv": (1,27),
                "germline_sv": (1,25),
                "germline_high_cov": (1,15),
                "somatic_tumor_only": (1,22),
                "somatic_fastpass": (1,23),
                "somatic_ensemble": (1,38),
                "somatic_sv": (1,14)
            },
            "longread_dnaseq": {
                "nanopore": (1,5),
                "revio": (1,14)
            },
            "methylseq": {
                "bismark": (1,18),
                "gembs": (1,20),
                "dragen": (1,18),
                "hybrid": (1,20)
            },
            "nanopore_covseq": {
                "default": (1,9),
                "basecalling": (1,12)
            },
            "rnaseq":{
                "stringtie": (1,21),
                "variants": (1,25),
                "cancer": (1,30)
            },
            "rnaseq_denovo_assembly": {
                "trinity": (1,24),
                "seq2fun": (1,5)
            },
            "rnaseq_light":{"no_protocol":(1,8)}
        }

        pipeline_data = valid_steps.get(pipeline, {})
        valid_range = pipeline_data.get(protocol, pipeline_data.get("default"))
        valid_start, valid_end = valid_range

        for part in step_range.split(','):
            part = part.strip()
            if '-' in part:
                try:
                    start, end = map(int, part.split('-', 1))
                except ValueError:
                    print(f"[ERROR] '{part}' not in the correct step range format.")
                    return False
                if start > end or start < valid_start or end > valid_end:
                    print(f"[ERROR] Range '{part}' is out of bounds.\nPlease enter a valid step range within these bounds: ({valid_start}-{valid_end}).")
                    return False
            else:
                try:
                    step = int(part)
                except ValueError:
                    print(f"[ERROR] '{part}' is not a number.")
                    return False
                if step < valid_start or step > valid_end:
                    print(f"[ERROR] Step '{step}' is out of bounds.\nPlease enter a valid step range within these bounds: ({valid_start}-{valid_end}).")
                    return False

        return True

#for testing
def main():
    print("\nWelcome to the GenPipes Wizard!")
    print("This tool will help you select the appropriate deployment method, pipeline, protocol, and/or construct the command to run GenPipes.")
    print ("Let’s begin!\n")
    start_json_file = "general_guide.json"
    start = Wizard(start_json_file)
    start.tree_traversal()

if __name__ == "__main__":
    main()