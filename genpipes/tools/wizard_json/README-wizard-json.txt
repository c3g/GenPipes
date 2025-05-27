# Organization of Wizard JSON Files (refer to GenPipes_Wizard.drawio)

## `general_guide.JSON`
This file contains the general questions that the wizard will ask the user to determine which guide they need help with. 
In cases where the user skips a guide, they will be asked to select their choice of deployment method/pipeline/protocol.

- Deployment guide  
- Pipeline guide  
- Protocol guide  
- Command guide  
  - Within this guide, the user can also follow the step guide if needed

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
      - protocol_help if pipeline requires a protocol
      - command_help if pipeline does not require a protocol

- protocol_help: Asks the user if they need help choosing a protocol
  - Yes → start of protocol_guide.json
  - No → protocol_selection

- protocol_selection: Presents protocol options based on the selected pipeline
  - Uses "choices_cases" with "when" clauses to filter protocol list
  - Selection leads to <protocol_name>_protocol_selected

- <protocol_name>_protocol_selected:
  - Stores protocol_name variable
  - Goes to command_help

- command_help: Asks if the user needs help constructing the command
  - Yes → start of command_guide.json
  - No → end

- end: Terminates the wizard with a final message
  - Suggests running "genpipes -h" or visiting ReadTheDocs for more support


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
