<!-- Improved compatibility of back to top link: See: https://github.com/othneildrew/Best-README-Template/pull/73 -->
<a id="readme-top"></a>
<!--
*** Thanks for checking out the Best-README-Template. If you have a suggestion
*** that would make this better, please fork the repo and create a pull request
*** or simply open an issue with the tag "enhancement".
*** Don't forget to give the project a star!
*** Thanks again! Now go create something AMAZING! :D
-->

<!-- PROJECT SHIELDS -->
<!--
*** I'm using markdown "reference style" links for readability.
*** Reference links are enclosed in brackets [ ] instead of parentheses ( ).
*** See the bottom of this document for the declaration of the reference variables
*** for contributors-url, forks-url, etc. This is an optional, concise syntax you may use.
*** https://www.markdownguide.org/basic-syntax/#reference-style-links
-->

<!-- PROJECT LOGO -->
<br />
<div align="center">
  </a>

<h3 align="center">GATK-For-Microbes Best Practice Workflow: BV-BRC Adaptation</h3>

  <p align="center">
    The GATK-For-Microbes Best Practice Workflow configured to work in the BV-BRC and ANL system via Snakemake.
    <br />
    <a href="https://github.com/alemanac/GATK_mutect2_shifted_module_ANL"><strong>Explore the docs »</strong></a>
    <br />
    <br />
    <a href="https://github.com/alemanac/GATK_mutect2_shifted_module_ANL/issues/new?labels=bug&template=bug-report---.md">Report Bug</a>
    &middot;
    <a href="https://github.com/alemanac/GATK_mutect2_shifted_module_ANL/issues/new?labels=enhancement&template=feature-request---.md">Request Feature</a>
  </p>
</div>

<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li><a href="#about-the-project">About The Project</a>
      <ul>
        <li><a href="#summary-and-key-features">Summary and Key Features</a></li>
        <li><a href="#differences-from-public-workflow">Differences From Public Workflow</a></li>
        <li><a href="#dag-and-summary">DAG and Summary</a></li>
        <li><a href="#directory-structure">Directory Structure</a></li>
        <li><a href="#format-of-output-files">Format of Output Files</a></li>
        <li><a href="#important-output-files">Important Output Files</a></li>
        <li><a href="#purpose-of-snakefiles">Purpose of Snakefiles</a></li>
      </ul>
    </li>
    <li><a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a>
      <ul>
        <li><a href="#quickstart-i-suggest-beginning-here">Quickstart</a></li>
        <li><a href="#snakemake-execution">Snakemake Execution</a></li>
        <li><a href="#advance-usage">Advanced Usage</a>
          <ul>
            <li><a href="#advance-config-parameters">Advanced Config Parameters</a></li>
            <li><a href="#advance-snakemake-execution">Advanced Snakemake Execution</a></li>
            <li><a href="#specifying-config-filepath">Specifying Config Filepath</a></li>
            <li><a href="#specfiying-config-parameters-through-the-command-line">Specifying Config Parameters Through the Command Line</a></li>
          </ul>
        </li>
      </ul>
    </li>
    <li><a href="#roadmap">Roadmap</a></li>
    <li><a href="#contributing">Contributing</a>
      <ul>
        <li><a href="#top-contributors">Top Contributors</a></li>
      </ul>
    </li>
    <li><a href="#faq">FAQ</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgments">Acknowledgments</a></li>
  </ol>
</details>


<!-- ABOUT THE PROJECT -->
## About The Project

### Summary and Key Features

This is a Snakemake module which follows the [GATK-For-Microbes Best Practice Workflow](https://gatk.broadinstitute.org/hc/en-us/articles/360060004292-Introducing-GATK-for-Microbes). This workflow implements the following novel features not seen in traditional GATK or variant-detection-based workflows:

* Shifting of reference FASTA to detect edge variants near the beginning and end of the reference FASTA.
* Resolution of "dangling ends" found when locally assembling reads.
* variant detection based on *both* local assembly and pileup data.
* Rescuing of true variants in heavily softclipped regions.

### Differences From Public Workflow

This workflow, configured for the BV-BRC and ANL system, differs from the public version in a couple of ways.

* No usage of containers. This workflow assumes you are executing Snakemake through one of the BV-BRC-ANL ubuntu containers. It is not recommended to launch a container from within a container, hence this difference.
* Shell commands are slightly different. Some of the containers used in the public version call commands differently than how we do. For example, `bwa` is called via `/usr/gitc/bwa` in a container used in a couple of the rules in the alignment process; however, we call `bwa` via just `bwa` in our ubuntu container. Differences like these had to be made.
* The usage of a precompiled GATK jre package. The public workflow does not need to a path towards a precompiled .jar for GATK, as it launches GATK via its containers. However, since our ubuntu container does not have GATK preinstalled nor aliased, this workflow requires you to manually download the precompiled GATK4 package. Fortunately, it is easy and instructions are described in the Installation section.

### DAG

#### DAG and Summary

The Snakemake module has the following DAG:

![DAG]

This may seem daunting, but essentially the DAG follows these steps:

1. Index the reference FASTA and shifted reference FASTA.
2. Align reads to both FASTAs.
3. Call variants against both reference FASTAs.
4. Merge vcf files.
5. Filter variants.

### Directory Structure

This workflow has the following directory structure:

* Note: this follows the Snakemake reproducibility guidelines.
* Another note: the results directory may not be used and/or created depending on how your `config.yaml` (or the config you specify; more details in the advance usage section) is configured. More on this in the Usage section.
* Another *another* note: the `resources` directory may not be used. it depends if you map your own reference FASTA file via a command line argument when executing Snakemake. This is described in the advance usage section.

```sh
.
├── dag.pdf
├── dag.png
├── README.md
├── resources
│   └── example_ref.fna
├── results
│   └── example_sample
│       └── example_rule
│           └── output_files.ext
└── workflow
    ├── config
    │   └── config.yaml
    ├── rules
    │   ├── alignment_shifted.smk
    │   ├── alignment.smk
    │   ├── mutect2.smk
    │   ├── ref_processing_shifted.smk
    │   └── ref_processing.smk
    └── snakefile
```

#### Format of Output Files

As described above, the output files of the workflow are organized as so, where `base_output_dir` can be configured in the `config.yaml` file:

* `base_output_dir/sample/rule_name/files`
* `base_output_dir/ref_base/ref_processing`
  * Every output file for each rule in ref_processing.smk is outputted here.
* `base_output_dir/ref_base/ref_processing_shifted`
  * Every output file for each rule in ref_processing_shifted.smk is outputted here.

#### Important Output Files

The following files are probably of high importance:

* `{main_dir}/{SRR}/remove_duplicate_rows/variants.vcf`: final .vcf file of the overall workflow.
  * Note: the unfiltered variants are still within the .vcf file, if interested.
* `{main_dir}/{SRR}/gen_duplicates_marked_reads/reads.bam`: the final state of the aligned reads used for variant calling via `Mutect2`.
* `{main_dir}/{ref_base}/ref_processing_shifted/{ref_base}.shifted{ref_ext}"`: the shifted reference FASTA file used for the shifted portion of the workflow.

#### Purpose of Snakefiles

Below is a summary of what each Snakefile posseses and does.

* `workflow/snakefile`: serves as the "master" snakefile which imports the other snakefiles.
* `workflow/rules/ref_processing.smk`: indexes the reference FASTA file. This snakefile contains the following rules:
  * `move_ref`
  * `gen_ref_dict`
  * `gen_ref_faidx`
  * `gen_ref_bwa_indexes`
* `workflow/rules/ref_processing_shifted.smk`: this is basically snakefile as `workflow/rules/ref_processing.smk`, but contains the `shift_fasta` rule and inherits the rules from ref_processing. By "inheriting the rules," I mean using the exact same shell but with different output and inputs. In other words, this snakefile indexes the shifted FASTA file. This snakefile contains the following rules:
  * `move_ref_shifted`
  * `shift_fasta`
  * `gen_ref_bwa_indexes_shifted`
* `workflow/rules/alignment.smk`: this aligns the reads of the SRR from the `samples` part of the config (described in the usage section) to your reference FASTA file. This snakefile contains the following rules:
  * `fetch_reads`
  * `gen_aligned_reads`
  * `gen_unmapped_bam`
  * `gen_merged_aligned_reads`
  * `gen_duplicates_marked_reads`
  * `gen_sorted_reads`
* `workflow/rules/alignment_shifted.smk`: this does the same thing as the snakefile above but aligns the reads to the shifted reference FASTA file instead. This snakemake contains the following rules:
  * `gen_aligned_reads_shifted`
  * `gen_merged_aligned_reads_shifted`
  * `gen_duplicates_marked_reads_shifted`
  * `gen_sorted_reads_shifted`
* `workflow/rules/mutect2.smk`: this workflow uses `Mutect2` for variant calling. It calls variants against the regular and shifted reference FASTA files. Then, it merges the variants produced from both vcf files (after de-offsetting the shifted variants) and filters them. It contains the following rules:
  * `gen_mutect2_vcfs`
  * `gen_mutect2_vcfs_shifted`
  * `lifted_over_and_combined_vcfs`
  * `merge_stats`
  * `gen_filtered_vcfs`
  * `remove_duplicate_rows`: this rule might not be useful anymore. I previously got duplicated rows from how I previously structured the workflow, and hence needed this rule – but it's been restructured, so perhaps this rule isn't warranted anymore.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- GETTING STARTED -->
## Getting Started

### Prerequisites

The module was tested to work correctly with ubuntu-075-12.sif and likely works with higher versions. The only necessary change is the version of java mapped mapped to `java`. GATK recommends java 8, but I believe I ran into an error with that version, so I've instead been using java 23. The installation section describes how I installed and mapped, via `sdkman`, `java` to java 23 in the ubuntu container.

Additionally, you'll need to download the GATK jar and set the absolute filepath to the config file used to run the Snakemake module (more on this in the usage section). If interested, the [Getting with GATK 4 docs](https://gatk.broadinstitute.org/hc/en-us/articles/360036194592-Getting-started-with-GATK4) describes how to set an alias to `GATK`. However, the Snakemake module will then need to be modified, as the shell commands operate with the pre-compiled .jar file set in the config file.

* `Java 23.0.1`
  * Which was installed via `SDKMAN!`. More details in the installation section.
* `GATK >= 4.6.0.0`
  * You only need the precompiled package. More details in the installation section.
* Everything else is already within the ubuntu container, which is the following:
  * `BWA`
  * `Samtools`
  * `fasterq-dump`
  * `Snakemake`

### Installation

1. Begin installing [SDKMAN!](https://sdkman.io) via the command

    ```sh
    curl -s "https://get.sdkman.io" | bash
    ```

2. Add the following to the bottom of your `.bashrc` (`SDKMAN!` should prompt you do to this):

    ```sh
    #THIS MUST BE AT THE END OF THE FILE FOR SDKMAN TO WORK!!!
    export SDKMAN_DIR="$HOME/.sdkman"
    [[ -s "$HOME/.sdkman/bin/sdkman-init.sh" ]] && source "$HOME/.sdkman/bin/   sdkman-init.sh"

    source "/home/ac.aleman/.sdkman/bin/sdkman-init.sh"
    ````

3. Using `SDKMAN!`, install `java 23.0.1` via

    ```sh
    sdk install 23.0.1-amzn
    ```

4. Initialize `SDKMAN!` (this must be done everytime you launch the container) via

   ```sh
   source "$HOME/.sdkman/bin/sdkman-init.sh"
   ```

5. Set `java` to `java 23.0.1` using `SDKMAN!` (this also must be done everytime you open the container) via

   ```sh
   sdk use java 23.0.1-amzn
   ```

6. Download the latest GATK package following step 4 of [Getting Started with GATK4](https://gatk.broadinstitute.org/hc/en-us/articles/360036194592-Getting-started-with-GATK4)

7. Make note of the full filepath to `gatk-package-[version]-local.jar`. Doing so will be important for the Usage section.

8. git clone the module

   ```sh
   git clone https://github.com/alemanac/Bacterial_GATK_SNPs_aleman.git
   ```

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- USAGE EXAMPLES -->
## Usage

### Quickstart (I suggest beginning here!)

#### Config Parameters

Quickly take a look at the `config.yaml` and make note of what you see.

To summarize, here the most important parts of `config.yaml`:

* `gatk_jar: "absolute filepath inside quotes"`: obtain the absolute filepath of the .jar file for GATK4, as noted in step 7, and put it here.
* `ref_fna: "absolute filepath"` or move and/or copy this into the resources section and just specify `"resources/ref_fna.ext"`: this is the reference FASTA used for variant-calling purposes. I recommend doing a full run of what I preconfigured, so do not change this just yet.
* This part spans multiple lines, so please reference the block code below. These are your samples to be analyzed for variants. I recommend just keeping this as so for the reasons above.
  * Note: unfortunately, I have yet to add support for a .bam file to be used directly or your own fastq files. This would be a trivial task with some Snakemake knowledge. If needed, please feel free to contact me, and I will gladly add support ASAP!

```yaml
samples:
- "SRR3722077"
- "SRR3722078"
```

Continuing on...

* `output_dir: "results"`: this is the directory where your output directories and files will be put. Think of it as the `{base_dir}` described in the Directory Structure section. It can be an absolute filepath or relative to the root of the workflow-directory structure. Right now, this is simply `results`, meaning all output directories and files will be outputed in the `results` directory as shown in the directory tree in the Directory Structure secton.
* `allele_fraction: "0.75"`: this is the minimum allele fraction needed for a variant to be filtered.
  * Note: variants can still fail to be filtered under other conditions.

Again, I recommend not changing anything except for the `gatk_jar` filepath on your first try – just to ensure if things work.

#### Snakemake Execution

Before running Snakemake, I recommend checking if the dag can be properly built given your config – in other words, if the wildcards passed from the `snakefile` can be passed to an output file from another rule, and if the wildcards there can be passed on to another output file from another rule, and so on.

You can do so via the following command:

`snakemake -np`

Caution: you must be within the root directory of the workflow. I recognize this may be an inconvenience (unfortunately a bit too late), but here is a quick hack to resolve this:

`CWD=(pwd); cd path_to_root_of_workflow; snakemake -np; cd $CWD`

This `cd`s you into the root of the workflow and `cd`s you back into your current working directory before you ran that part of the command.

If no errors are met, you can run the snakemake workflow via the following:

`snakemake --cores num` where num is the number of cores you want to use. Snakemake can run jobs in parallel, so keep that in mind when specifying this value.

Or, if you want to resume the previous hack:

`CWD=(pwd); cd path_to_root_of_workflow; snakemake --cores num; cd $CWD`

### Advance Usage

#### Advance Config Parameters

Here are some more advance config features.

* `java_args: "-Xms4000m -Xmx8000m"`: this denotes how much memory you allocate to the heap of the JRE. `-Xms` denotes the initial amount and `-Xmx` denotes the final amount.
* See the block below, but this is metadata related to the reads. This shouldn't be relevant to the BV-BRC system, I think, but I provided them as arguments in the config for public use – and kept them here in the ANL version of this workflow.

```yaml
read_group_library: "NA"
read_group_platform: "NA"
read_group_platform_unit: "NA"
read_group_run_date: "NA"
read_group_sequencing_center: "NA"
```

#### Advance Snakemake Execution

##### Specifying config filepath

If you're creating your own pipeline utilizing this Snakemake workflow, you might want to specify your own configfile. You can do via

`snakemake --configfile  path/to/config.yml --other_args...`

Note: you will still need to cc do the root of the directory as described above and copied below:

`CWD=(pwd); cd path_to_root_of_workflow; snakemake --configfile  path/to/config.yml --other_args...; cd $CWD`

Another note: By not specifying this config path, Snakemake will assume you are using `workflow/config/config.yaml`.

##### Specfiying Config Parameters Through The Command Line

You may find it inconvenient to change the config parameters manually or via a script. Fortunately, you can change them via the CLI like so

`snakemake --config allele_frac="1.5"`

I honestly haven't used this feature that much, so that above example may not work, but I hope it gives the general idea. The official example of this feature on the Snakemake docs can be found [here](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html).

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- ROADMAP -->
## Roadmap

See the [open issues](https://github.com/alemanac/GATK_mutect2_shifted_module_ANL/issues) for a full list of proposed features (and known issues).

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- CONTRIBUTING -->
## Contributing

I used a template for this README, and below is what it recommends – but I'm flexible!

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

<p align="right">(<a href="#readme-top">back to top</a>)</p>

### Top contributors

<a href="https://github.com/alemanac/GATK_mutect2_shifted_module_ANL/graphs/contributors">
  <img src="https://contrib.rocks/image?repo=alemanac/GATK_mutect2_shifted_module_ANL" alt="contrib.rocks image" />
</a>

## FAQ

todo...


<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- CONTACT -->
## Contact

Slack me @ Tony Aleman!

Project Link: [https://github.com/alemanac/GATK_mutect2_shifted_module_ANL](https://github.com/alemanac/GATK_mutect2_shifted_module_ANL)

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- ACKNOWLEDGMENTS -->
## Acknowledgments

* [Allan Dickerman](https://github.com/AllanDickerman), for his guidance.
* Rebecca Wattam, also, for her guidance.
* [GATK-For-Microbes, by the Broad Institute](https://github.com/broadinstitute/GATK-for-Microbes).
* The open source community – for helping with my debugging!

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[DAG]: dag.png