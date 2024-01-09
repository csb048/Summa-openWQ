# Install Summa-OpenWQ

summa-openwq link https://github.com/csb048/Summa-openWQ/tree/file_splitting_revised_install

openwq link https://github.com/csb048/openwq/tree/develop_revised_install

### Installation instuctions derived from...
https://simcity.usask.ca/trac/ticket/3207

## Docker

Ensure that Docker is installed on your machine. LINK TO DOCKER

## Installation

### Obtaining Scource Code
* Clone Summa-Openwq -> https://github.com/csb048/Summa-openWQ.git
    > Use the file_splitting_revised_install branch until updated, remember to remove this line
    
    `git clone https://github.com/csb048/Summa-openWQ.git`

    This will add a folder containg the Summa souce code into the current working dirrectory

* Clone OpenWQ -> https://github.com/csb048/openwq.git
    * cd into `Summa-openWQ/build/source/openwq`
        
        `cd Summa-openWQ/build/source/openwq`
    * Clone OpenWq into the `openwq` folder
        > Use the develop_revised_install branch until updated, remember to remove this line
    
        `git clone https://github.com/csb048/openwq.git` 

### The Docker Enviroment

In `Summa-openWQ/utils/containers/summa_openwq/` are the files needed for creating and launching the docker container.

*  Open launch_docker_container.sh in a text editor
    * Change the `PROJECT_DIR` variable to the dirrectory containing the vesion of `Summa-OpenWQ` that was cloned (not `OpenWQ`).
    * Change the `DATA_DIR` variable to a dirrectory contiaing the data you want to use for running the simulations
    
    Now these two dirrectories can be accses and modified from within the container.

* Run the `build_docker_container.sh` script to build the container
    `./build_docker_container.sh`

* Now run the `launch_docker_container.sh` script **AFTER** modifying it to match your enviroment
    `./launch_docker_container.sh`

* Attach to the Docker Enviroment
    * Use `docker ps` to get infomation about the container.  Use `docker attach <CONTAINER ID>` to attach via shell.
    * Attach to the container in VS Code using the Docker extentions... 
        * Docker: https://marketplace.visualstudio.com/items?itemName=ms-azuretools.vscode-docker
        * Dev Containers: https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers

* Once you ace attached to the container cd into `/code`, there you will see the two dirrectories we mapped ealier, `Summa-OpenWQ/` and `data/`.

### Compiling Summa-OpenWQ
* We are going to make some modifications to opwnwq's `CMakeLists.txt` found at `/code/Summa-openWQ/build/source/openwq/openwq/CMakeLists.txt`. Open `CMakeLists.txt` with a text editor.
    * Ensure that `COMPILE_TARGET` is set to `summa_openwq`.
    * Set `COMPILE_MODE` to `debug` or `fast` depending on your needs.
        * `degug` is good for developing and making changes to the softwear
        * `fast` is for if you want to just use the program.          
* cd into `/code/Summa-OpenWQ/build/`
* Create a new directory called `cmake_build`
* cd into `cmake_build/`
* run `../call_cmake.sh`
    * If you have any issues double check the paths in the call_cmake.sh file.

* run `make -j 2`
    * This should create a summa_openwq executable in `/code/case_studies`

If you come across Issues with compiling try deleting everything in `cmake_build/` and try running `../call_cmake.sh` again, Make sure you are in the correct directory and the directory in `call_cmake.sh` leads to `../source/openwq/openwq/.

### Setup Synthetic Tests
Before any simulation is ran we have to mofify the configuration files to match docker enviroment. For any simulation you have to modify the Summa and OpenWQ config files.  As an example We will be configuring a synthetic test. Clone the `synthetic_tests` repo -> LINK HERE into `data/`.  **USE ABSOLUTE PATHS** 
* Summa
    * Modify the fileManager `summa_fileManager_OpenWQ_systheticTests_BGQ.txt` located in `/code/data/synthetic_tests/13_batch_oxygenBODcycle/summa/summa/SUMMA`, Ensure the settingsPath, forcingPath, and output feilds lead to the correct location.
        * settingsPath: `/code/data/synthetic_tests/13_batch_oxygenBODcycle/summa/summa/SUMMA/`
        * forcingPath:  `/code/data/synthetic_tests/13_batch_oxygenBODcycle/summa/summa/forcing_data/`
        * outputPath:   `/code/data/synthetic_tests/13_batch_oxygenBODcycle/summa/summa/output/`

        When changing the output feild, ensure that the directory it leads to already exists.  It is very important that each of these directories have a "/" at the end.
* OpenWQ
    * inside `/code/data/synthetic_tests/13_batch_oxygenBODcycle/summa/` there are several `.json` files, open `openWQ_master.json`.  Ensure that each Json file that contains a FILEPATH entry points to the correct location.
        * `"CONFIG_FILEPATH": "/code/data/synthetic_tests/13_batch_oxygenBODcycle/summa/openWQ_config_13.json"`
        * `"MODULE_CONFIG_FILEPATH": "/code/data/synthetic_tests/13_batch_oxygenBODcycle/summa/openWQ_MODULE_BGC_cycling_13.json"`
        * `"MODULE_CONFIG_FILEPATH": "/code/data/synthetic_tests/13_batch_oxygenBODcycle/summa/openWQ_MODULE_TE_batchSimulations.json"	`
        
        You can change `RESULTS_FOLDERPATH` to somewhere specific while editing `openWQ_master.json`.
        * `"RESULTS_FOLDERPATH": "/code/data/synthetic_tests/13_batch_oxygenBODcycle/summa/summa/Output_OpenWQ"`

This example only applies for Synthetic test #13.  Modifiying config for other datasets follows the same general steps.

### Running Summa-OpenWQ

* Make an enviroment variable containing the absolute path of the master json file being used.
    `export master_json=PATH/TO/MASTER/JSON`

* `./summa_openwq_debug -g 1 1 -m path/to/fileManager.txt`

# Common Issues
* When Running Docker commands (Including the scripts for building and launching docker containers) you may get an error stating **"ERROR: permissions denied ..."**.  To fix this run the problematic command with **sudo**,  `sudo <rest of command>`. or add youself to the docker group on your local machine.

*  When building the container some of the apt-gets will fail for whatever reason with error code 100.  to fix this, cry about it ;_;  (no solution yet).  Ive only came accros this issue on my personal laptop.  I've tried used the updated Dockerfile that kyle made Late Dec 2023 but it still tries to install the problematic dependancies.



