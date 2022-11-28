
<details>
    <summary> <summary>
        <style>
            H2{color:RoyalBlue !important;}
            H3{color:DarkOrange !important;}
            H4{color:#16E92C !important;}
            p{color:White !important;}
        </style>
</details>

# SEG-Filter Algorithm
---
## **Table of Content**
---
1. Introduction 

2. Scripts or Programs
    * SEG Algorithm
    * `FASTA_conv.py`
    * `seg_scr.sh`
    * `main.py`
    * `visual.py`
    * assets folder 
    * `statistics.py `

3. Getting Started
    * Dependencies
    * Installation
    * Formating

4. Script Execution and Results
    * `seg_scr.sh`
    * `main.py`
    * `visual.py` and assests folder
    * `statistics.py`
<br /><br />

----
## **Scripts or Program** 
----

#### SEG Algorithm

* We will be using the 21.10 version of the ncbi-seg version of the SEG algorithm, this program will identify and mask low-complexity regions within amino acid sequences. 

<br /> 

#### `seg_scr.sh`

* This bash script functions:
    * run the SEG Program on all of the FASTA files inside a folder

    * output each FASTA files into their own individual .txt files 

* purpose: 
    * save user from having to run the SEG program on each FASTA files manually.

#### `main.py`

* This python script is the SEG-Filter Algorithm.

#### `visual.py` and assets folder

* This python script functions:
    * Produce the cluster graphs 

    * The assets folder is necessary for the stylaztion of the graphs that the `visual.py` needs to output nice graphs. 

* Purpose:
    * Produce cluster graphs to visualize any significant trends within the dataset. 

#### `statistics.py`
* This python script functions:
    * Calculate _Conditional Probabilities_

    * Calculate _Multinomial Distributon_ for each parsed sequence using the known [probability distribution of all amino acids.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7127678/pdf/main.pdf)


* Purpose:
    * Analyze the likelihood of each sequence occuring randomly.

<br />

----
## **Getting Started**
----

### **Dependencies**
1. It is required that the user has updated their Python to atleast version 3.7.0 to run any of the provided python scripts correctly. 

2. python libraries required
    * collections 
    * argparse
    * pandas 
    * math
    * random
    * PIL 
    * assets.circlify 
    * matplotlib.pyplot  


### **Installation** 

#### **SEG Algorithm**
1.  Download the 21.10 version of the [ncbi-seg Algorithm](http://manpages.ubuntu.com/manpages/impish/man1/ncbi-seg.1.html).
 
2. Configuration 
    * Window User
        * configure the program normally in your linux enviroment.
    * macOS User
        * __Option 1__: Install Window support software through Bootcamp Assitant using this [guide](https://support.apple.com/boot-camp). This method will cost the user at least 64 gb of storage space. 

        * __Option 2__: If possible, it is recocommended the user configure the SEG Program on a separate linux manchine and remotely connect to the server from your local machine. 

#### **Scripts**
1. The rest of the python and bash scripts can be found in the [Script Folder](https://ucdavis.app.box.com/folder/135515480079) of the murraylab box.
    * `seg_scr.sh`
    * `main.py`
    * `visual.py`
    * assests folder 
    * `sequence.py`

#### **Formatting SEG Input Files**
1. The formatting of the description line for FASTA files should follow this specific rule:

**Example**
```
Formatting rule:
>"sequence name" [family:"family it comes from"]

Example: 
“>PRNP [family:prp]
```

* Deviation from the above format will result in the main.py script erroring out when running it on the wrongly formatted SEG output .txt files. 
    * Note:
        * The users can change the keywords the `main.py` script searches for in the description line when parsing through the seg output files. This can be done by manipulating a few lines of codes within the `main.py` script. 
        
        * Any manipulation of the key words within the `main.py` script should also be done in the `statistics.py` script as well. 

2. The descript line contains two important components:

    1. sequence name

    2. the family the sequence comes from. 

3. The SEG-Filter Algorithm only search for two things in the description line: 
    1. the string attached to > 
        * if the string is multiple words long then connect the words together with the underscore function.

            * ex) Major histocompatibility antigens needs to become Major_histocompatibility_antigens

    2. the string that cames after "family:"

**Example**
``` 
 >PRNP [family:prp]

 1. the string attached to > 

>PRNP 

the string would be "PRNP"

2. the string that cames after "family:"

[family:prp]

the string would be "prp"
``` 

<br />

---
## **Script Execution and Results**
----
* Run the scripts based the following order. 

### **1) seg_scr.sh**

* Run the `seg.scr.sh`  on the same directory level as the SEG program and the input folder with all the FASTA files. 

![example](https://github.com/Truc-T-Le/SEG-Algorithm-Project/blob/main/Example1.jpg)


**Code**
```
Generic: 
./seg_scr.sh "input folder" "output folder"

Example:
./seg_scr.sh input_folder seg_output
```

**Results**
* Running the `seg_scr.sh` script will give you an output folder (Using the example above, the folder name is called "seg_output") that contains all the output text files from the SEG Algorithm. 

![example](https://github.com/Truc-T-Le/SEG-Algorithm-Project/blob/20e90972c9873fcfafe0357d7487c420344ba779/Result1.png)


**Notes for macOS user**
* macOS creates hidden /.ds files within the "seg_output" folder once you run the ./seg_scr.sh bash script on a host server and copy the folder onto your local server, these files will cause the `main.py` script to error out.

* To avoid the next step from crashing, run the following two command lines:

```
# To find if there was anything ./ds file
for f in `find .`; do echo `file -I "$f"`; done

# To remove the files
rm -v **/.ds
```

<br />

### **2) `main.py`**

* Run the `main.py` on the same directory level as the SEG output folder. 

![example](https://github.com/Truc-T-Le/SEG-Algorithm-Project/blob/18223f062bc46f153c2c2d9007d1217f163b595a/Example2.png)

***Code***
```
Generic:
Python3 main.py "input folder"

Example:
Python3 main.py seg_output
```

**Results**

* Running the `main.py` script will output a folder called <span style="color:red">**output**</span>.

![Part 1](https://github.com/Truc-T-Le/SEG-Algorithm-Project/blob/0000b2cc1ef0bd0ace29496b82b868cabbff179c/Result2_1.png) 

<br />

* The **output** folder contains three other folders:
    * <span style="color:red">**cluster**</span> folder
        * contains multiple csv files based on the top n amino acids component within each sequence. 
    * <span style="color:red">**graphs**</span> folder
        * Contains 3 addtional folders:
            * <span style="color:red">by_family</span> folder
            * <span style="color:red">by_id</span> folder 
            * <span style="color:red">by_overall</span> folder
    * <span style="color:red">**sequences**</span> folder 
        * Contains 3 additional folders:
            * <span style="color:red">by_family</span> folder
            * <span style="color:red">by_id</span> folder 
            * <span style="color:red">by_overall</span> folder

![Part 2](https://github.com/Truc-T-Le/SEG-Algorithm-Project/blob/0000b2cc1ef0bd0ace29496b82b868cabbff179c/Result2_2.png)

    
### **3) `visual.py` and Assest folder**

* The `visual.py` and the <span style="color:red">**assest**</span> folder should be located in the same directory as the <span style="color:red">**output**</span> folder before you run the `visual.py` script.
    * **Note!!!**: 
        *  To avoid the `visual.py` python script from erroring out, you need to have the <span style="color:red">**assest**</span> folder in the same directory!!!


![Example](https://github.com/Truc-T-Le/SEG-Algorithm-Project/blob/8b03c5586cfee664af3bd658a899a1d887178c37/Example3.png)

**Code**

```
python3 visual.py output
```

**Result**

* Running the `visual.py` script will output 3 clusters graphs:
    * gc-output.jpg
    * gcs-output.jpg
    * gs-output.jpg

![Example](https://github.com/Truc-T-Le/SEG-Algorithm-Project/blob/d45d4273022b219b27a3b26124c4c10664ad4232/Result3_1.png)


### **4) `statistics.py`**

* Run the `statistics.py` script in the same directory as the <span style="color:red">**output**</span> folder.

![Example](https://github.com/Truc-T-Le/SEG-Algorithm-Project/blob/5e2ee5cb29ffb3231b2c32d6392b4938c369c529/Example4.png)

**Code**

```
python3 statistics.py output
```

**Result**

* Running the `statistics.py` script will output a folder called <span style="color:red">**probability**</span> in the <span style="color:red">**output**</span> folder. 

![Example](https://github.com/Truc-T-Le/SEG-Algorithm-Project/blob/399f946a89688ae73faf44f501c68440007b3492/Result4.png)

* The <span style="color:red">**probability**</span> folder contains:
    * <span style="color:red">**by_family**</span> folder
        * <span style="color:red">**order_doesn't_matter**</span> folder
        
        * <span style="color:red">**order_does_matter**</span> folder
            
    * <span style="color:red">**by_id**</span> folder
        * <span style="color:red">**order_doesn't_matter**</span> folder

        * <span style="color:red">**order_does_matter**</span> folder

    * <span style="color:red">**by_overall**</span> folder
        * <span style="color:red">**order_doesn't_matter**</span> folder

        * <span style="color:red">**order_does_matter**</span> folder

    * `multinomial.csv` 

![Example](https://github.com/Truc-T-Le/SEG-Algorithm-Project/blob/399f946a89688ae73faf44f501c68440007b3492/Result4_2.png)

![Example](https://github.com/Truc-T-Le/SEG-Algorithm-Project/blob/399f946a89688ae73faf44f501c68440007b3492/Result4_3.png)

**Notes**
* <span style="color:red">**order_doesn't_matter**</span> folder: 
    * Neglect the positional ordering of the amino acids within the pairings.
        * Ex) AT and TA are consider the same pairing of amino acids,

    * Recorded counts of occurance for each possible pairing of amino acids.

    * Calculated the conditional probability for each possible pairs of amino acids. 


* <span style="color:red">**order_does_matter**</span> folder: 
    * Take into consideration the positional ordering of the amino acids within the pairings.
        * EX) AT and TA are consider unique pairings of amino acids.

    * Recorded counts of occurance foe each unique possible pairing of amino acids.
    
    * Calculated the conditional probability for each unique possible pairing of amino acids.

* The calculated conditonal probabilities can be useful for constructing a _Hidden Markov's Model_ to calculate the likeliness of each sequence occuring. 


## **Notes**

### The main.py script 
* the script can be modified to accomodate the users needs.
    * HOWEVER
        * If the user modified the cluster function to output more than the top 2 amino acids, which is the current script default function, they would need to adjust the visual.py script accordingly as well. 

### The statistics.py script 
* This script can only calculate the conditional probabilities of two amino acids pairings (2 events).

* As of July 6th, 2021, the author does not have the require mathematic knowledge to construct an algorithm that calculates the conditional probabilities with multiple conditions. (She probably can, but is currently refusing to accept the harsh reality.)
