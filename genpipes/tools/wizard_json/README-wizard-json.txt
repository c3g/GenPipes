# Organization of Wizard JSON Files (refer to GenPipes_Wizard.drawio)

## `general_guide.JSON`
This file contains the general questions that the wizard will ask the user to determine which guide they need help with. 
In cases where the user skips a guide, they will be asked to select their choice of deployment method/pipeline/protocol.

- Deployment guide  
- Pipeline guide  
- Protocol guide  
- Command guide  
  - Within this guide, the user can also follow the step guide if needed

Legend of the node names and their function:

- start_general_guide: Asks the user if they need help deploying GenPipes
  - Yes → goes to the start of deployment_guide.json
  - No → goes to pipeline_help

- pipeline_help: Asks if the user needs help selecting the appropriate pipeline
  - Yes → goes to the start of pipeline_guide.json
  - No → goes to pipeline_selection

- pipeline_selection: Lets the user select a pipeline
  - Selection leads to <pipeline_name>_pipeline_selected

- <pipeline_name>_pipeline_selected:
  - Stores pipeline_name variable
  - Goes to:
      - protocol_help if pipeline requires a protocol
      - command_help if pipeline does not require a protocol

- protocol_help: Asks the user if they need help choosing a protocol
  - Yes → goes to the start of protocol_guide.json
  - No → goes to protocol_selection

- protocol_selection: Presents protocol options based on the selected pipeline
  - Uses "choices_cases" with "when" clauses to filter protocol list
  - Selection leads to <protocol_name>_protocol_selected

- <protocol_name>_protocol_selected:
  - Stores protocol_name variable
  - Goes to command_help

- command_help: Asks if the user needs help constructing the command
  - Yes → goes to the start of command_guide.json
  - No → goes to end

- end: Terminates the wizard with a final message
  - Suggests running "genpipes -h" or visiting ReadTheDocs for more support


## `deployment_guide.JSON`
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


## `pipeline_guide.JSON`
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
