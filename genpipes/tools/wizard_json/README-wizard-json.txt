# Organization of Wizard JSON Files (refer to GenPipes_Wizard.drawio)

##List of variables in nodes
- id
- type
- question
- options
- label
- next
- external: name of json file to enter
- entryPoint: id of node to enter
- node: id of node to go to next
- message
- variable: name of variable 
- value: information stored in variable
- choices_cases
- when ... equals ... pipeline_name
- choices 
- cases
- prompt: text to ask user for input

##List of types 
- confirm: yes/no questions
- selection: multi-select questions
- set_variable: store variable
- message: output message for user
- switch
- input: input from user 

##List of set_variable
From general_guide.JSON:
- pipeline_name: ampliconseq, chipseq, covseq, dnaseq, longread_dnaseq, methylseq, nanopore_covseq, rnaseq, 
                  rnaseq_denovo_assembly, rnaseq_light
- protocol_name: chipseq, atacseq, germline_snv, germline_sv, germline_high_cov, somatic_tumor_only, somatic_fastpass,
                 somatic_ensemble, somatic_sv, nanopore, revio, bismark, gembs, dragen, hybrid, stringtie, variants, cancer, 
                 trinity, seq2fun

From command_guide.JSON:
- r_command: -r {1_TODO_IN_PYTHON}
- j_command: -j slurm, -j pbs, -j batch
- scheduler_server_name: beluga, cedar, narval, abacus, batch 
- server_in: GENPIPES_INIS/common_ini/{scheduler_server_name}.ini
- path_custom_ini: 2_TODO_IN_PYTHON
- c_command: -c $GENPIPES_INIS/{pipeline_name}/{pipeline_name}.base.ini
             -c $GENPIPES_INIS/{pipeline_name}/{pipeline_name}.base.ini $GENPIPES_INIS/common_ini/{scheduler_server_name}.ini
             -c $GENPIPES_INIS/{pipeline_name}/{pipeline_name}.base.ini $GENPIPES_INIS/{pipeline_name}/{path_custom_ini}.ini
             -c $GENPIPES_INIS/{pipeline_name}/{pipeline_name}.base.ini $GENPIPES_INIS/common_ini/{scheduler_server_name}.ini $GENPIPES_INIS/{pipeline_name}/{path_custom_ini}.ini
- o_command: {empty placeholder if user skips}, -o {directory_name}
- d_command: {empty placeholder if user skips}, -d {design_file_name}.{pipeline_name}.txt
- p_command: {empty placeholder if user skips}, -p {pair_file_name}.dnaseq.txt
- s_command: {empty placeholder if user skips}, -s {step_range}
- g_command: -g {3_TODO_IN_PYTHON}
- final_command: genpipes {pipeline_name} -t {protocol_name} {c_command} {r_command} {d_command} {p_command} {j_command} {s_command} {o_command} {g_command}

From step_guide.JSON:
- step_range: 1-6, 8 ; 1-17, 19-23 ; 1-18, 20-24 ; 1-14, 17-18 ; 1-16, 19-20 ; 1-18, 20-21 ; 1-19, 21-24 ; 1-3, 5 ; 1-6, 8

##TODO_IN_PYTHON
- 1st TODO_IN_PYTHON: insert the readset filename (depending if input contains .txt)
- 2nd TODO_IN_PYTHON: insert the path to custom ini (depending if input contains .ini)
- 3rd TODO_IN_PYTHON: insert the genpipes filename (depending if input contains .sh)

------------------------------------------------------------------------------------------------------------------------------------
## `general_guide.JSON`
This file contains the general questions that the wizard will ask the user to determine which guide they need help with. 
In cases where the user skips a guide, they will be asked to select their choice of deployment method/pipeline/protocol.

- Deployment guide  
- Pipeline guide  
- Protocol guide  
- Command guide (within this guide, the user can also follow the step guide if needed)

**Legend of the node names and their function:**
- start_general_guide: Asks the user if they need help deploying GenPipes
  - Yes → start of deployment_guide.json
  - No → pipeline_help

- pipeline_help: Asks if the user needs help selecting the appropriate pipeline
  - Yes → start of pipeline_guide.json
  - No → pipeline_selection

- pipeline_selection: Lets the user select a pipeline
  - Selection leads to <pipeline_name>_pipeline_selected

- <pipeline_name>_pipeline_selected:
  - Stores pipeline_name variable
  - Goes to:
      - pipeline_choice_message_next_protocol if pipeline requires a protocol
      - pipeline_choice_message_next_command if pipeline does not require a protocol

- pipeline_choice_message_next_protocol: Display pipeline choice
  - goes to protocol_help

- pipeline_choice_message_next_command: Display pipeline choice
  - goes to command_help

- protocol_help: Asks the user if they need help choosing a protocol
  - Yes → start of protocol_guide.json
  - No → protocol_selection

- protocol_selection: Presents protocol options based on the selected pipeline
  - Uses "choices_cases" with "when" clauses to filter protocol list
  - Selection leads to <protocol_name>_protocol_selected

- <protocol_name>_protocol_selected:
  - Stores protocol_name variable
  - goes to protocol_choice_message

- protocol_choice_message: Display protocol choice
  - goes to command_help
  
- command_help: Asks if the user needs help constructing the command
  - Yes → start of command_guide.json
  - No → end

- end: Terminates the wizard with a final message
  - Suggests running "genpipes -h" or visiting ReadTheDocs for more support

------------------------------------------------------------------------------------------------------------------------------------
## `deployment_guide.JSON`
This file contains the questions that help the user determine the deployment method they want to use to deploy GenPipes

**Deployment method options:**
- `DRAC infrastucture`
- `cloud`
- `container`
- `locally`

**Legend of the node names and their function:**
- start_deployment_guide: Asks user if they have access to DRAC
  - Yes → drac_deployment_selected
  - No → cloud_question

- cloud_question: Asks user if they want to use cloud infrastructure
  - Yes → cloud_deployment_selected
  - No → container_question

- container_question: Asks user if they want to use a local cluster and need access the CVMFS
  - Yes → container_deployment_selected
  - No → local_deployment_selected

- drac_deployment_selected: suggests to follow link for detailed deployment steps
  - goes to pipeline_help

- cloud_deployment_selected: suggests to follow link for detailed deployment steps
  - goes to pipeline_help

- container_deployment_selected: suggests to follow link for detailed deployment steps
  - goes to pipeline_help

- local_deployment_selected: suggests to follow link for detailed deployment steps
  - goes to pipeline_help

------------------------------------------------------------------------------------------------------------------------------------

## `pipeline_guide.JSON`
This file contains the questions that help the user determine the appropriate pipeline based on their dataset and analysis goals.  

**Pipeline options:**
- `ampliconseq`
- `chipseq`
- `covseq`
- `dnaseq`
- `longread_dnaseq`
- `methylseq`
- `nanopore_covseq`
- `rnaseq`
- `rnaseq_denovo_assembly`
- `rnaseq_light`

**Legend of the node names and their functions:**
- start_pipeline_guide: Asks whether the dataset originates from DNA samples or from RNA/COVID samples
  - Yes → chip_dna_longread_methyl_seq_question  
  - No → rnaseq_question

- chip_dna_longread_methyl_seq_question: Distinguishes between long-read methylation/long-read DNA and methylation/ChIP/DNA sequencing pipelines
  - Yes → longread_methylseq_question  
  - No → methylseq_question

- longread_methylseq_question: Determines whether the user wants to analyze long-read methylation vs long-read DNA 
  - Yes → longread_methylseq_pipeline_selected  
  - No → longread_dnaseq_pipeline_selected

- methylseq_question: Determines whether the user wants to analyze DNA methylation or proceed to ChIP/DNA sequencing 
  - Yes → methylseq_pipeline_selected  
  - No → chip_dna_seq_question

- chip_dna_seq_question: Chooses between ChIP-seq and DNA-seq pipelines
  - Yes → chipseq_pipeline_selected  
  - No → dnaseq_pipeline_selected

- rnaseq_question: Distinguishes between RNA-seq pipelines and nanopore/cov/amplicon pipelines 
  - Yes → 3_rnaseq_question  
  - No → covid_question

- 3_rnaseq_question: Determines whether the RNA dataset is from a non-model organism, selecting between the 3 RNA-seq options
  - Yes → rnaseq_denovo_assembly_pipeline_selected  
  - No → 2_rnaseq_question

- 2_rnaseq_question: Chooses between the 2 remaining RNA-seq pipelines 
  - Yes → rnaseq_pipeline_selected  
  - No → rnaseq_light_pipeline_selected

- covid_question: Distinguishes between nanopore/cov and amplicon sequencing pipelines  
  - Yes → nanopore_cov_seq_question  
  - No → ampliconseq_pipeline_selected

- nanopore_cov_seq_question: Chooses between nanopore and cov sequencing pipelines  
  - Yes → nanopore_covseq_pipeline_selected  
  - No → covseq_pipeline_selected

Note: All <pipeline_name>_pipeline_selected nodes jump to the appropriate entry point node in general_guide.json.

------------------------------------------------------------------------------------------------------------------------------------

## `protocol_guide.JSON`
This file contains the questions used to determine the appropriate protocol based on the dataset and analysis goals.
Note: if ampliconseq/nanopore_covseq/covseq/rnaseq_light then skip question asking user if they need pipeline guide.

**Protocol options:**

- `chipseq` → `chipseq`, `atacseq`  
- `dnaseq` → `germline_snv`, `germline_sv`, `germline_high_cov`, `somatic_tumor_only`, `somatic_fastpass`, 
`somatic_ensemble`, `somatic_sv`  
- `longread_dnaseq` → `nanopore`, `revio`  
- `methylseq` → `bismark`, `gembs`, `dragen`, `hybrid`  
- `rnaseq` → `stringtie`, `variants`, `cancer`  
- `rnaseq_denovo_assembly` → `trinity`, `seq2fun`

**Legend of the node names and their functions:**
- start_protocol_guide: Entry point that branches based on selected pipeline  
  - chipseq → chipseq_protocol_help  
  - dnaseq → dnaseq_protocol_help  
  - longread_dnaseq → longread_dnaseq_protocol_help  
  - methylseq → methylseq_protocol_help  
  - rnaseq → rnaseq_protocol_help  
  - rnaseq_denovo_assembly → rnaseq_denovo_assembly_protocol_help

- chipseq_protocol_help: Asks whether the user wants to analyze DNA-protein interactions  
  - Yes → chipseq_protocol_selected  
  - No → atacseq_protocol_selected

- dnaseq_protocol_help: Asks whether the analysis focuses on germline DNA variants  
  - Yes → germline_question  
  - No → somatic_question

- germline_question: Asks whether dataset has high coverage for detecting low-frequency variation  
  - Yes →  germline_snv_protocol_selected  
  - No → germline_sv_snv_question

- germline_sv_snv_question: Asks whether to focus on structural variants (SV) or single nucleotide variants (SNV)
  - SV → germline_sv_protocol_selected  
  - SNV → germline_snv_protocol_selected  

- somatic_question: Asks if dataset has matched tumor and normal samples
  - Yes → somatic_ensemble_protocol_selected  
  - No → tumor_only_question

- tumor_only_question: Asks if dataset only has tumor
  - Yes → somatic_tumor_only_protocol_selected  
  - No → somatic_sv_question

- somatic_sv_question: Asks if need to perform somatic structural variant detection  
  - Yes → somatic_sv_protocol_selected  
  - No → somatic_fastpass_protocol_selected

- longread_dnaseq_protocol_help: Asks if data is generated on Nanopore instrument
  - Yes → nanopore_protocol_selected  
  - No → revio_protocol_selected

- methylseq_protocol_help: Asks if user has access to Illumina Dragen  
  - Yes → dragen_hybrid_protocol_help  
  - No → bismark_gembs_protocol_help

- dragen_hybrid_protocol_help: Lets user choose between Dragen and Hybrid protocols  
  - dragen → dragen_protocol_selected  
  - hybrid → hybrid_protocol_selected

- bismark_gembs_protocol_help: Lets user choose between Bismark and GemBS protocols  
  - bismark → bismark_protocol_selected  
  - gembs → gembs_protocol_selected

- rnaseq_protocol_help: Asks if the dataset comes from cancer samples  
  - Yes → cancer_protocol_selected  
  - No → variant_question

- variant_question: Asks if variant detection is required  
  - Yes → variants_protocol_selected  
  - No → stringtie_protocol_selected

- rnaseq_denovo_assembly_protocol_help: Asks if organism lacks a reference genome  
  - Yes → seq2fun_protocol_selected  
  - No → trinity_protocol_selected

Note: All <protocol_name>_protocol_selected nodes jump to the appropriate entry point node in general_guide.json.

------------------------------------------------------------------------------------------------------------------------------------

## `command_guide.JSON`
This file contains the questions that the wizard will ask the user to construct the appropriate command based on their pipeline, 
protocol, readset file, job scheduler, design/pair file, directory, steps, etc.

- Pipeline  
- Protocol  
- Readset file  
- Job scheduler  
- Design/pair file  
- Output directory  
- Steps, etc.

**Legend of the node names and their functions:**

- start_command_guide: Asks the user to enter the readset file name
  - store_readset_filename

- store_readset_filename: Stores readset filename as part of r_command
  - slurm_job_scheduler_question

- slurm_job_scheduler_question: Asks if the user wants to use SLURM
  - Yes → slurm_j_command
  - No → general_job_scheduler_question

- slurm_j_command: Sets SLURM as the job scheduler (-j slurm)
  - slurm_server_question

- slurm_server_question: Asks user to select SLURM server
  - Beluga → beluga_server_ini
  - Cedar → cedar_server_ini
  - Narval → narval_server_ini

- beluga_server_ini: Stores 'beluga' as the scheduler_server_name
  - store_name_for_server_ini

- cedar_server_ini: Stores 'cedar' as the scheduler_server_name
  - store_name_for_server_ini

- narval_server_ini: Stores 'narval' as the scheduler_server_name
  - store_name_for_server_ini

- general_job_scheduler_question: Asks user to select a general job scheduler
  - Abacus → abacus_j_command
  - Batch → batch_j_command

- abacus_j_command: Sets PBS for Abacus cluster (-j pbs)
  -  abacus_server

- abacus_server: Sets scheduler_server_name to 'abacus'
  - custom_ini_question

- batch_j_command: Sets Batch scheduler (-j batch)
  - batch_server_ini

- batch_server_ini: Sets scheduler_server_name to 'batch'
  - store_name_for_server_ini

- store_name_for_server_ini: Constructs path to server-specific ini file
  - custom_ini_question

- custom_ini_question: Branches based on pipeline_name
  - pipeline where certain protcols have custom ini:
    - dnaseq → custom_ini_check_protocol
  - pipeline with custom ini:
     - nanopore_covseq → o_command_construction_with_custom_choice
  - other pipelines : c_command_construction_no_custom_choice

- custom_ini_check_protocol: Checks protocol_name for pipelines that support custom ini's
  - protocols_with_custom:
    - germline_sv → o_command_construction_with_custom_choice
    - somatic_fastpass  → o_command_construction_with_custom_choice
    - somatic_ensemble  → o_command_construction_with_custom_choice
    - somatic_sv  → o_command_construction_with_custom_choice
  - other protocols → c_command_construction_no_custom_choice

- path_custom_ini_question: Asks user for path to custom ini
  - store_path_custom_ini

- store_path_custom_ini: Stores and formats path to custom ini
  - o_command_construction_with_custom_choice

- c_command_construction_no_custom_choice: Sets default -c command based on scheduler_server_name
  - abacus → c_command_construction_abacus_no_custom
  - other scheduler servers → c_command_construction_no_custom

- c_command_construction_abacus_no_custom: Uses only pipeline base ini (no server ini)
  - initialize_o_command

- c_command_construction_no_custom: Includes both pipeline and server ini paths
  - initialize_o_command

- c_command_construction_abacus_with_custom: Sets -c using pipeline base ini and custom ini (abacus)
  - initialize_o_command

- c_command_construction_with_custom: Sets -c using base ini, server ini, and custom ini
  - initialize_o_command

- initialize_o_command: Initializes o_command
  - o_command_question

- o_command_question: Asks if user wants to change working directory
  - Yes → input_directory_name
  - No → pipeline_check_for_file_question

- input_directory_name: Asks user for new directory name
  - store_directory_name

- store_directory_name: Sets -o with user directory
  - initialize_d_command

- initialize_d_command: Initializes d_command
  - initialize_p_command

- initialize_p_command: Initializes p_command
  - pipeline_check_for_file_question

- pipeline_check_for_file_question: Based on pipeline_name, determine next input requirement
  - Pipeline has design file:
    - ampliconseq
    - chipseq
    - methylseq
    - rnaseq (protocol: stringtie)
    - rnaseq_denovo_assembly
    - rnaseq_light
  - Pipeline has pair file:
    - dnaseq (somatic_tumor_only protocol has no pair file )
  - Pipeline has no design/pair file:
    - covseq
    - longread_dnaseq
    - nanopore_covseq

- design_file_question: Asks if the user has a design file
  - Yes → input_design_file_name
  - No → start_step_guide in step_guide.json

- input_design_file_name: Prompts for design file name
  - store_design_file_name

- store_design_file_name: Formats and sets -d using design file
  - step_question

- check_protocol_pair_file: Checks protocol_name to determine if pair file is required
  - tumor_only → pair_file_question
  - somatic_fastpass, somatic_ensemble, somatic_sv → input_pair_file_name
  - others → step_question

- pair_file_question: Asks if user has a pair file
  - Yes → input_pair_file_name
  - No → step_question

- input_pair_file_name: Prompts for pair file name
  - store_pair_file_name

- store_pair_file_name: Formats and sets -d with pair file for dnaseq
  - initialize_s_command

- initialize_s_command: Initializes s_command
  - step_question

- step_question: Asks if user wants to run all pipeline steps
  - Yes → g_command_question
  - No → input_step_range

- input_step_range: Prompts for specific steps to run
  - store_step_range

- store_step_range: Sets the -s argument with the provided step range
  - g_command_question

- g_command_question: Asks user to enter file name to run GenPipes
  - store_g_file_name

- store_g_file_name: Sets -g with g command file name 
  - command_construction

- command_construction: Construst full command to run GenPipes
  - Replaces placeholders: "genpipes {pipeline_name} -t {protocol_name} {c_command} {r_command} {d_command} {p_command} {j_command} {s_command} {o_command} {g_command}
  - goes to end 

- end: Displays command

------------------------------------------------------------------------------------------------------------------------------------

## `step_guide.JSON`
This file defines the step range to be added to the command based on the pipeline and protocol, 
specifically when the user does not have a design file and wants to run GenPipes. Steps involving the design file must be skipped.

**Step-based options:**

- `ampliconseq`  
- `chipseq` → `chipseq`, `atacseq`  
- `methylseq` → `bismark/hybrid`, `gembs/dragen`  
- `rnaseq` → `stringtie`  
- `rnaseq_denovo_assembly` → `trinity`, `seq2fun`  
- `rnaseq_light`

**Legend of the node names and their functions:**
- start_step_guide: Entry point that branches based on the selected pipeline
  - ampliconseq → ampliconseq_step_command  
  - chipseq → chipseq_step_help  
  - methylseq → methylseq_step_help  
  - rnaseq → rnaseq_step_help  
  - rnaseq_denovo_assembly → rnaseq_denovo_assembly_step_help  
  - rnaseq_light → rnaseq_light_step_command

- ampliconseq_step_command: Sets step range to "1-6, 8"
  - goes to initialize_g_command in command_guide.json

- chipseq_step_help: Entry point that branches based on selected protocol
  - chipseq → chipseq_step_command  
  - atacseq → atacseq_step_command

- chipseq_step_command: Sets step range to "1-17, 19-23"
  - goes to initialize_g_command in command_guide.json

- atacseq_step_command: Sets step range to "1-18, 20-24"
  - goes to initialize_g_command in command_guide.json

- methylseq_step_help: Entry point that branches based on selected protocol
  - bismark → bismark_hybrid_step_command  
  - hybrid → bismark_hybrid_step_command  
  - dragen → dragen_gembs_step_command  
  - gembs → dragen_gembs_step_command

- bismark_hybrid_step_command: Sets step range to "1-14, 17-18"
  - goes to initialize_g_command in command_guide.json

- dragen_gembs_step_command: Sets step range to "1-16, 19-20"
  - goes to initialize_g_command in command_guide.json

- rnaseq_step_command: Sets step range to "1-18, 20-21"
  - goes to initialize_g_command in command_guide.json

- rnaseq_denovo_assembly_step_help: Entry point that branches based on selected protocol
  - trinity → trinity_step_command  
  - seq2fun → seq2fun_step_command

- trinity_step_command: Sets step range to "1-19, 21-24"
  - goes to initialize_g_command in command_guide.json

- seq2fun_step_command: Sets step range to "1-3, 5"
  - goes to initialize_g_command in command_guide.json

- rnaseq_light_step_command: Sets step range to "1-6, 8"
  - goes to initialize_g_command in command_guide.json