---
title: "Data & Setup"
number-sections: false
---

<!-- 
Note for Training Developers:
We provide instructions for commonly-used software as commented sections below.
Uncomment the sections relevant for your materials, and add additional instructions where needed (e.g. specific packages used).
Note that we use tabsets to provide instructions for all three major operating systems.
-->

::: {.callout-tip level=2}
## Workshop Attendees

If you are attending one of our workshops, we will provide a training environment with all of the required software and data.  
If you want to setup your own computer to run the analysis demonstrated on this course, you can follow the instructions below.
:::

## Data

The data used in these materials is not yet publicly available. 
We will add a link to the data in due time.

## Software

### R and RStudio

::: {.panel-tabset group="os"}
#### Windows

Download and install all these using default options:

- [R](https://cran.r-project.org/bin/windows/base/release.html)
- [RTools](https://cran.r-project.org/bin/windows/Rtools/)
- [RStudio](https://www.rstudio.com/products/rstudio/download/#download)

#### Mac OS

Download and install all these using default options:

- [R](https://cran.r-project.org/bin/macosx/)
- [RStudio](https://www.rstudio.com/products/rstudio/download/#download)

#### Linux

- Go to the [R installation](https://cran.r-project.org/bin/linux/) folder and look at the instructions for your distribution.
- Download the [RStudio](https://www.rstudio.com/products/rstudio/download/#download) installer for your distribution and install it using your package manager.
:::

#### R Packages

Open RStudio. 
In the R console, run the following commands to install all the necessary packages: 

```r
install.packages("BiocManager")
BiocManager::install(c("dada2", "phyloseq", "Biostrings", "ggplot2", "reshape2", "readxl", "tidyverse"))
```


### Linux {#sec-install-linux}

:::{.panel-tabset}
#### Fresh Installation

The recommendation for bioinformatic analysis is to have a dedicated computer running a Linux distribution. 
The kind of distribution you choose is not critical, but we recommend **Ubuntu** if you are unsure. 

You can follow the [installation tutorial on the Ubuntu webpage](https://ubuntu.com/tutorials/install-ubuntu-desktop#1-overview). 

:::{.callout-warning}
Installing Ubuntu on the computer will remove any other operating system you had previously installed, and can lead to data loss. 
:::

#### Windows WSL

The **Windows Subsystem for Linux (WSL2)** runs a compiled version of Ubuntu natively on Windows. 

There are detailed instructions on how to install WSL on the [Microsoft documentation page](https://learn.microsoft.com/en-us/windows/wsl/install). 
But briefly:

- Click the Windows key and search for  _Windows PowerShell_, right-click on the app and choose **Run as administrator**. 
- Answer "Yes" when it asks if you want the App to make changes on your computer. 
- A terminal will open; run the command: `wsl --install`.  
  Progress bars will show while installing "Virtual Machine Platform", "Windows Subsystem for Linux" and finally "Ubuntu" (this process can take a long time).
    - **Note:** it has happened to us in the past that the terminal freezes at the step of installing "Ubuntu". If it is frozen for ~1h at that step, press <kbd>Ctrl + C</kdb> and hopefully you will get a message saying "Ubuntu installed successfully".
- After installation completes, restart your computer.
- After restart, a terminal window will open asking you to create a username and password.  
  If it doesn't, click the Windows key and search for _Ubuntu_, click on the App and it should open a new terminal. 
  - You can use the same username and password that you have on Windows, or a different one - it's your choice. Spaces and other special characters are not allowed for your Ubuntu username.
  - **Note:** when you type your password nothing seems to be happening as the cursor doesn't move. However, the terminal is recording your password as you type. You will be asked to type the new password again to confirm it, so you can always try again if you get it wrong the first time.

You should now have access to a Ubuntu Linux terminal. 
This behaves very much like a regular Ubuntu server. 

##### Configuring WSL2

After installation, it is useful to **create shortcuts to your files on Windows**. 
Your main `C:\` drive is located in `/mnt/c/` and other drives will be equally available based on their letter. 
To create shortcuts to commonly-used directories you use _symbolic links_. 
Here are some commands to automatically create shortcuts to your Windows "Documents",  "Desktop" and "Downloads" folders (copy/paste these commands on the terminal):

```bash
ln -s $(wslpath $(powershell.exe '[environment]::getfolderpath("MyDocuments")' | tr -d '\r')) ~/Documents
ln -s $(wslpath $(powershell.exe '[environment]::getfolderpath("Desktop")' | tr -d '\r')) ~/Desktop
ln -s $(wslpath $(powershell.exe '[environment]::getfolderpath("UserProfile")' | tr -d '\r'))/Downloads ~/Downloads
```

You may also want to **configure the Windows terminal to automatically open _WSL2_** (instead of the default Windows Command Prompt or Powershell):

- Search for and open the "<i class="fa-solid fa-terminal"></i> Terminal" application.
- Click on the down arrow <i class="fa-solid fa-chevron-down"></i> in the toolbar.
- Click on "<i class="fa-solid fa-gear"></i> Settings".
- Under "Default Profile" select "<i class="fa-brands fa-linux"></i> Ubuntu".


#### Virtual Machine

Another way to run Linux within Windows (or macOS) is to install a Virtual Machine.
However, this is mostly suitable for practicing and **not suitable for real data analysis**.

Details for installing Ubuntu on VirtualBox is given on [this page](https://ubuntu.com/tutorials/how-to-run-ubuntu-desktop-on-a-virtual-machine-using-virtualbox#1-overview).
Make sure to do these things, while you are setting it up:

- In Step 2 "Create a user profile": make sure to tick the Guest Additions option.
- In Step 2 "Define the Virtual Machine’s resources": 
  - Assign at least 4 CPUs and 16000MB of RAM. At the very minimum you need 2 CPUs to run an Ubuntu VM.
  - Set at least 100GB as disk size, more if you have it available (note, this will not take 100GB of space on your computer, but it will allow using up to a maximum of that value, which is useful as we are working with sequencing data).

Once the installation completes, login to the Ubuntu Virtual machine, open a terminal and do the following: 

- Run `su` command.
- Enter your user password. Your terminal should change to start with `root@`
- Type the command: `usermod -a -G sudo YOUR-USERNAME-HERE`.
- Close the terminal and restart the virtual machine. 

These commands will add your newly created user to the "sudo" (admin) group. 
:::


After making a fresh install of Ubuntu (using any of the methods above), open a terminal and run the following commands to update your system and install some essential packages: 

```bash
sudo apt update && sudo apt upgrade -y && sudo apt autoremove -y
sudo apt install -y git
sudo apt install -y default-jre
```


### _Conda_ 

We recommend using the _Conda_ package manager to install your software. 
In particular, the newest implementation called _Mamba_. 

To install _Mamba_, run the following commands from the terminal: 

```bash
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh -b -p $HOME/miniforge3
rm Miniforge3-$(uname)-$(uname -m).sh
$HOME/miniforge3/bin/mamba init
```

Restart your terminal (or open a new one) and confirm that your shell now starts with the word `(base)`.
Then run the following commands: 

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set remote_read_timeout_secs 1000
```


### Bioinformatics Software

We can install all the software with `mamba`:

```bash
mamba create -n metagen

mamba install -n metagen fastqc multiqc cutadapt trimmomatic bowtie2 samtools metaphlan mash SPAdes bbmap flash maxbin2 prokka gtdbtk abricate checkm-genome
```

From now on, you can use these packages, by activating the software environment using `mamba activate metagen`.


### Databases

Some of the programs used require us to download public databases in addition to their installation. 
To follow these instructions, make sure you have activated the software environment first: `mamba activate metagen`.

#### CheckM (275MB)

The [CheckM documentation](https://github.com/Ecogenomics/CheckM/wiki/Installation#required-reference-data) gives the link to its database file. 

We will download this databases to a directory called `checkmdb`, but you can change this if you prefer to save it elsewhere. 

```bash
mkdir checkmdb
```

Download and decompress the file:

```bash
wget -O checkm_db.tar.gz https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
tar -xzvf checkm_db.tar.gz -C checkmdb
rm checkm_db.tar.gz
```

After downloading, we run the following command to configure CheckM:

```bash
checkm data setRoot $(pwd)/checkmdb/
```

#### GTDB-Tk (40GB)

This program offers a script to automatically download the database:

```bash
download-db.sh
```
