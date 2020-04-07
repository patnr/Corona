<!--Corona-->
<!------------>

This repo holds research on the Corona epidemic. 
<!--[Read the Docs](https://coronastudies.readthedocs.io/en/latest/)-->


![Screenshot from epcalc.py](./screenshot_epcalc.png)

Get started
================================================
Works on Linux/Windows/Mac.

1. **Prerequisite**: Python>=3.7.  
   If you're not a python expert:  
   1a. Install Python via [Anaconda](https://www.anaconda.com/download).  
   1b. Use the [Anaconda terminal](https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html#starting-conda) to run the commands below.  
   1c. (Optional) [Create & activate a new Python environment](https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html#managing-environments).  
   If the installation (below) fails, try doing step 1c first.

2. **Install**:  
   `$ git clone https://github.com/patricknraanes/Corona.git`  
   `$ pip install -e Corona`  
   NB: don't download via your browser (instead of git). This zip won't install.

3. **Test:**  
   `$ cd Corona/src/corona`  
   `$ python epcalc.py`  
   You should get output like the figure above.  

PS: personally, I prefer to use ipython:  
   `$ ipython`  
   `In [1]: run epcalc.py`  
