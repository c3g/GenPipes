#!/usr/bin/env python3

import json
import os
import sys
import logging

import questionary 
from jinja2 import Environment 

# Configure logging at the module level
wiz_log = logging.getLogger("wizard")
wiz_log.setLevel(logging.INFO)
# Define a handler
console_handler = logging.StreamHandler(sys.stdout)
console_handler.setLevel(logging.INFO)
# Define a formatter
console_formatter = logging.Formatter("%(message)s")
wiz_log.addHandler(console_handler)
console_handler.setFormatter(console_formatter)

def load_guide (file_path):
    """
    Load wizard JSON files by filename: {general, deployment, pipeline, protocol, command, step}_guide.json
    """
    full_file_path = os.path.join(os.path.dirname(__file__), "wizard_json", file_path)
    with open(full_file_path) as file:
        return json.load(file)
    
class Wizard:

    def __init__ (self, start_file):
        #store variables
        self.variables = {}

        #Jinja2 environment for evaluating placeholders {{...}}
        self.env = Environment()

        #load initial guide
        self.current_guide = load_guide(start_file)

        #track current json file being used 
        self.current_file = start_file

        #set first node to start traversal
        self.current_node_id = self.current_guide["_meta"]["entry_point"]

        #keep track of node id's visited 
        self.history = []

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
        current = self.find_node(self.current_node_id)
        if current.get("type") in ("input", "confirm", "selection"):
            self.history.append((self.current_file, self.current_node_id))

        # when next node is in the current JSON file
        if isinstance(next_node, str):
            self.current_node_id = next_node
        else:
            #when next node is in another JSON file
            self.current_file = next_node["external"]
            self.current_guide = load_guide(self.current_file)
            self.current_node_id = next_node["entryPoint"]

    def go_back(self):
        """
        Go back to previous question
        """
        while self.history:
            self.current_file, self.current_node_id = self.history.pop()
            self.current_guide = load_guide(self.current_file)
            node = self.find_node(self.current_node_id)

            # Clear value of variable if type is input or set_variable 
            if node.get("type") in ("input", "set_variable") and "variable" in node:
                self.variables.pop(node["variable"], None)

            if node.get("type") != "message":
                return True

        wiz_log.error("ERROR! You cannot go back further, you're at the first question.")
        return False
        
    def exit_text (self, prompt):
        """
        Prompt user to input text. Exit wizard if user cancels (Ctrl +C)
        """
        answer = questionary.text(f"{prompt} (or type 'back' to go back):").ask()
        if answer is None:
            wiz_log.info("Exiting GenPipes wizard.")
            sys.exit(0)
        return answer.strip()
    
    def exit_confirm(self, prompt):
        """
        Prompt user to answer Y/n question. Exit wizard if user cancels (Ctrl +C)
        """
        answer = questionary.select(prompt, choices=["Yes", "No", "Back"]).ask()
        if answer is None:
            wiz_log.info("\nExiting GenPipes wizard.")
            sys.exit(0)
        return answer
    
    def exit_select (self, prompt, choices):
        """
        Prompt user to select from list. Exit wizard if user cancels (Ctrl +C)
        """
        extended_choices = choices + ["Back"]
        answer = questionary.select(prompt, choices=extended_choices).ask()
        if answer is None:
            wiz_log.info("\nExiting GenPipes wizard.")
            sys.exit(0)
        return answer
    
    def add_space (self):
        """
        Add empty line between answer and the next question 
        """
        print()
    
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
                while True:
                    answer = self.exit_confirm(self.apply_variables(node["question"]))
                    if answer == "Back":
                        if self.go_back():
                            self.add_space()
                            break   
                        else:
                            continue
                    chosen = answer
                    next_info = next(opt["next"] for opt in node["options"] if opt["label"] == chosen)
                    self.goto(next_info)
                    self.add_space()
                    break

            #Selection: single-select question from list 
            elif node_type == "selection":
                #choices
                if "choices" in node:
                    labels = [current_choice["label"] for current_choice in node["choices"]]
                    while True:
                        choice = self.exit_select(self.apply_variables(node["question"]), labels)
                        if choice == "Back":
                            if self.go_back():
                                self.add_space()
                                break
                            else: 
                                continue   
                        next_node = next(current_choice for current_choice in node["choices"] if current_choice["label"] == choice)
                        self.goto(next_node["node"])
                        self.add_space()
                        break

                #choices_cases
                else:  
                    for case_block in node["choices_cases"]:
                        variable_name, value = next(iter(case_block["when"]["equals"].items()))
                        if self.variables.get(variable_name) == value:
                            labels = [current_choice["label"] for current_choice in case_block["choices"]]
                            while True:
                                choice = self.exit_select(self.apply_variables(node["question"]), labels)
                                if choice == "Back":
                                    if self.go_back(): 
                                        self.add_space()
                                        break
                                    else: 
                                        continue
                                next_node = next(current_choice for current_choice in case_block["choices"] if current_choice["label"] == choice)
                                self.goto(next_node["node"])
                                self.add_space()
                                break
                            break

            #Set_variable: set and store variable 
            elif node_type == "set_variable":
                variable = node["variable"]
                raw_value = node["value"]

                if variable in ("r_command", "path_custom_ini", "g_command", "d_command", "p_command", "genome_config_ini", "o_command"):
                    self.fix_filenames()


                updated_value = self.apply_variables(raw_value)

                #Clean up extra spaces in final command
                if variable == "final_command":
                    updated_value = " ".join(updated_value.split())

                self.variables[variable] = updated_value
                self.goto(node["next"])
                #self.add_space()

            #Message: output message for user, if no next node then end of wizard
            elif node_type == "message":
                message = self.apply_variables(node["message"])
                wiz_log.info(message)
                next_node = node.get("next")
                if not next_node:
                    break
                self.goto(next_node)
                self.add_space()

            #Switch: determine the next node based on cases (e.g pipeline/protocol name)
            elif node_type == "switch":
                variable = node["variable"]
                value = self.variables.get(variable)
                cases = node["cases"]
                if value in cases:
                    self.goto(cases[value]["node"])
                else:
                    wiz_log.error(f"ERROR! No matching case for {variable} ='{value}' at node {node['id']}")
                    sys.exit(1)

            #Input: Prompt the user for input and store it as a variable 
            elif node_type == "input":
                variable = node["variable"]
                while True:
                    input = self.exit_text(self.apply_variables(node["prompt"]))
                    if input.lower() == "back":
                        if self.go_back():
                            self.add_space() 
                            break
                        else: 
                            continue

                    self.variables[variable] = input.strip()

                    #ensure that user doesn't leave input questions empty
                    if variable == "raw_readset_filename":
                        if not input.strip():
                            self.add_space()
                            wiz_log.error("ERROR! You must enter a readset file name. Please try again.")
                            self.add_space()
                            continue
                    if variable == "design_file_name":
                        if not input.strip():
                            self.add_space()
                            wiz_log.error("ERROR! You must enter a design file name. Please try again.")
                            self.add_space()
                            continue
                    if variable == "pair_file_name":
                        if not input.strip():
                            self.add_space()
                            wiz_log.error("ERROR! You must enter a pair file name. Please try again.")
                            self.add_space()
                            continue
                    if variable == "raw_path_custom_ini":
                        if not input.strip():
                            self.add_space()
                            wiz_log.error("ERROR! You must enter a path to the custom ini name. Please try again.")
                            self.add_space()
                            continue
                    if variable == "directory_path":
                        if not input.strip():
                            self.add_space()
                            wiz_log.error("ERROR! You must enter a full path name to the directory. Please try again.")
                            self.add_space()
                            continue
                    if variable == "g_filename":
                        if not input.strip():
                            self.add_space()
                            wiz_log.error("ERROR! You must enter a file name. Please try again.")
                            self.add_space()
                            continue
                    if variable == "reference_genome_ini_filename":
                        if not input.strip():
                            self.add_space()
                            wiz_log.error("ERROR! You must enter a file name. Please try again.")
                            self.add_space()
                            continue

                    if variable == "step_range":
                        if not self.valid_step_range():
                            continue

                    self.goto(node["next"])
                    self.add_space()
                    break

            else:
                wiz_log.error(f"ERROR! Unknown node type: {node_type} in {self.current_file}")
                self.add_space()
                sys.exit(1)
    
    def fix_filenames(self):
        """
        Handle cases where user includes/doesn't include .txt/ini/sh to their input
        """
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
        if not self.variables.get("directory_path"):
            o_command = ""
        self.variables["o_command"] = o_command

        reference_genome_ini_filename = self.variables.get("reference_genome_ini_filename", "").strip()
        if reference_genome_ini_filename and not reference_genome_ini_filename.endswith(".ini"):
            reference_genome_ini_filename += ".ini"
        self.variables["reference_genome_ini_filename"] = reference_genome_ini_filename

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
        if not valid_range:
            self.add_space()
            wiz_log.error(f"ERROR! Please enter a valid step range.")
            self.add_space()
            return False
        valid_start, valid_end = valid_range

        for part in step_range.split(','):
            part = part.strip()
            if '-' in part:
                try:
                    start, end = map(int, part.split('-', 1))
                except ValueError:
                    self.add_space()
                    wiz_log.error(f"ERROR!'{part}' not in the correct step range format.")
                    self.add_space()
                    return False
                if start > end or start < valid_start or end > valid_end:
                    self.add_space()
                    wiz_log.error(f"ERROR! Range '{part}' is out of bounds.\nPlease enter a valid step range within these bounds: ({valid_start}-{valid_end}).")
                    self.add_space()
                    return False
            else:
                try:
                    step = int(part)
                except ValueError:
                    self.add_space()
                    wiz_log.error(f"ERROR! '{part}' is not a number.")
                    self.add_space()
                    return False
                if step < valid_start or step > valid_end:
                    self.add_space()
                    wiz_log.error(f"ERROR! Step '{step}' is out of bounds.\nPlease enter a valid step range within these bounds: ({valid_start}-{valid_end}).")
                    self.add_space()
                    return False

        return True

#for testing
def main():
    wiz_log.info(r"""
    __        __   _                            _          _   _                                               .-+: 
    \ \      / /__| | ___ ___  _ __ ___   ___  | |_ ___   | |_| |__   ___                                   .-+*--*- 
     \ \ /\ / / _ \ |/ __/ _ \| '_ ` _ \ / _ \ | __/ _ \  | __| '_ \ / _ \                                 :+#=.   =*:
      \ V  V /  __/ | (_| (_) | | | | | |  __/ | || (_) | | |_| | | |  __/                               :=#+.   ... -*:
      _\_/\_/ \___|_|\___\___/|_| |_| |_|\___| _\__\___/ __\__|_| |_|\___|    _ _                     .-#*.      :*=. =*: 
     / ___| ___ _ __ |  _ \(_)_ __   ___  ___  \ \      / (_)______ _ _ __ __| | |                 .-#*.       .=#+=#*=**:
    | |  _ / _ \ '_ \| |_) | | '_ \ / _ \/ __|  \ \ /\ / /| |_  / _` | '__/ _` | |              .:##:        .=#=     :*%+. 
    | |_| |  __/ | | |  __/| | |_) |  __/\__ \   \ V  V / | |/ / (_| | | | (_| |_|              ++.          =%- 
     \____|\___|_| |_|_|   |_| .__/ \___||___/    \_/\_/  |_/___\__,_|_|  \__,_(_)              .#-         .+=   
                             |_|                                                                -*:            -#.
                                                                                                :#*%@%=:  :+%@@*.-#:
    The GenPipes Wizard will help you select an appropriate deployment method,                   *+.            .++
    pipeline, protocol, and/or construct the command to run GenPipes.                    :*%%%%*=**:  ..:.....    .+#-+#%%#*:
                                                                                    .#%+:...    .*=                -#-     ..:+%#:  
    Press Ctrl+C at any time to exit the wizard.                                  #+.          .+##=.          .=#%%-          .=#.    
                                                                                %=              .:-=*%%%%%#+=-:.              -%.
    Let's begin!                                                                  .-*#=.                                    .=##=
                                                                                    .:=+*#**+=-:.                .:-=+**#*+=:. 
                                                                                            ..:---====+=++++=+===----:..
    """)                                                            
    start_json_file = "general_guide.json"
    start = Wizard(start_json_file)
    start.tree_traversal()

if __name__ == "__main__":
    main()
