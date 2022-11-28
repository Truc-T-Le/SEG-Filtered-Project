
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
    * Produce the cluster graphs and two CSVs that summarizes the graphs.

    * The assets folder is necessary for the stylaztion of the graphs that the `visual.py` needs to output nice graphs. 

* Purpose:
    * Produce cluster graphs to visualize any significant trends within the dataset. 

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
    * csv
    * shutil


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
1. The formatting of the description line for FASTA files should follow the conventional FASTA formatting [guidance](https://en.wikipedia.org/wiki/FASTA_format#Description_line) :

**Example**
```
Formatting rule:
>"information of the protein"

Examples: 
“>PRNP"
">sp|P31483|TIA1_HUMAN Cytotoxic granule associated RNA binding protein TIA1 OS=Homo sapiens OX=9606 GN=TIA1 PE=1 SV=3"
```

* Deviation from the above format will result in the main.py script erroring out. 


2. The SEG-Filter Algorithm only search for two things in the description line: 
    1. The ">" at the very start of the line. 
        
    2. The first string attached to the ">" symbol. 
        
3. The main.py script will take the first string between the ">" and the first white space or until the end of the description line.
        * The string is use to identify the protein that the LCR sub-sequence being analyze originated from.

         
**Example**

* ex) ">sp|P31483|TIA1_HUMAN Cytotoxic granule associated RNA binding protein TIA1 OS=Homo sapiens OX=9606 GN=TIA1 PE=1 SV=3" &#8594; "sp|P31483|TIA1_HUMAN"
        
* ex) “>PRNP" &#8594; “PRNP"


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

* Run the `main.py` on the same directory level as the folder contaning the SEG .txt files. 

![example](https://github.com/Truc-T-Le/SEG-Filtered-Project/blob/main/seg_result_1.png)

***Code***
```
Generic:
Python3 main.py "input folder"

Example:
Python3 main.py inputs
```

**Results**

* Running the `main.py` script will output a folder called <span style="color:red">**output**</span>.

![Part 1](https://github.com/Truc-T-Le/SEG-Filtered-Project/blob/main/seg_result_2.png) 

<br />

* The **output** folder contains three other folders:
    * <span style="color:red">**cluster**</span> folder
        * contains multiple csv files based on the top n amino acids component within each LCR sub-sequence. 
        * This example looks at the top 2 most abundant amino acids within each LCR sub-sequence.
        ![example](https://github.com/Truc-T-Le/SEG-Filtered-Project/blob/main/cluster.png) 
    * <span style="color:red">**graphs**</span> folder
        * Contains 2 addtional folders:
            * <span style="color:red">by_id</span> folder 
            * <span style="color:red">by_overall</span> folder
    * <span style="color:red">**sequences**</span> folder 
        * Contains 2 additional folders:
            * <span style="color:red">by_id</span> folder 
            * <span style="color:red">by_overall</span> folder

![Part 2](https://github.com/Truc-T-Le/SEG-Filtered-Project/blob/main/seg_result_3.png)

    
### **3) `visual.py` and Assest folder**

* The `visual.py` and the <span style="color:red">**assest**</span> folder should be located in the same directory as the <span style="color:red">**output**</span> folder before you run the `visual.py` script.
    * **Note!!!**: 
        *  To avoid the `visual.py` python script from erroring out, you need to have the <span style="color:red">**assest**</span> folder in the same directory as both the 'visual.py' and <span style="color:red">**output**</span> folder!!!


![Example](https://github.com/Truc-T-Le/SEG-Filtered-Project/blob/main/seg_result_2.png)

**Code**

```
python3 visual.py output
```

**Result**

* Running the `visual.py` script will output 2 clusters graphs and 2 CSV files:
    * gcs-output.jpg
    * gs-output.jpg
    * Networks_seq.csv
    * Networks.csv

![Example](https://github.com/Truc-T-Le/SEG-Filtered-Project/blob/main/visual_result.png)

## **Notes**

### The main.py script 
* The script can be modified to accomodate the users needs.
    * HOWEVER
        * If the user modified the cluster function to output more than the top 2 amino acids, which is the current script default function, they would need to adjust the visual.py script accordingly as well. 
        * User can modify this feature by modifying the first line within the ***def get_cluster_filenames*** in the 'main.py' script (line 518). 
        * Replace n in num=n with the desire number of amino acids within the pairing. 
        
        ```
        def get_cluster_filenames(args, common_letters_list):
        
            fmt_letters_list = list(map(lambda x : format_common_letters(x, num=n), common_letters_list))
            clusters = []
            for ll in fmt_letters_list:
                top_letters = list(map(lambda x: x[0], ll))
                clusters.append("".join(sorted("".join(top_letters))))

            f_clusters = []
            for cluster_name in clusters:
                f_cluster = os.path.join(args.out, "clusters", cluster_name + ".csv")
                f_clusters.append(f_cluster)

            return f_clusters

        ```
        
        * Within the 'visual.py' script, user needs to modify the pair object accordingly to the n most ambundant amino acids specify within the 'main.py' script. The current default is two amino acids, which is specify as pair[0] and pair[1] within the 'visual.py' script.

### The visual.py script
* The size of the graph outputted by the 'visual.py' script can be change by modifyng the following line of code within the def draw() within the Network class:
        ```
        width = height = int(len(self.group_nodes) ** a) * b
        ```
        * Recommended sizing: 
            * For Runs containing >50 LCR sub-sequence: a = 1.05, b = 120
            * For Runs containing 50 < x < 100: a = 1.15, b = 120
            * For Runs containing 100 < x < 150: a = 1.20, b = 125
            * For Runs containing  x < 500: a = 1.40, b = 145 
                                           
* Font size for the text in the graphs can be modify by changing the number within the following line:
          ```
          font = ImageFont.truetype("assets/Roboto-Light.ttf", 30)      
          ```
                                              
        
