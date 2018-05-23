# virtual-crossmatch

[VICTOR: VIrtual CrossmaTch for mOleculaR typing data](http://www.transplanttoolbox.org/victor)

This tool computes virtual crossmatch for molecular HLA typing data of donors (Genotype List Strings and NMDP Multiple Allele Codes) and gives the probabilities of positivity of crossmatch for conflicting antigens.

To use the user interface of the tool visit http://www.transplanttoolbox.org/victor

```
Steps to installing web app locally:
1. Install Python 3
2. Install Django, Django Rest Swagger and Django Rest Framework python modules
    - pip3 install django
    - More detailed instructions here: https://docs.djangoproject.com/en/1.11/topics/install/

3. Install Requests python module
    - pip3 install requests  
4. Clone this GitHub repository to a local directory (git clone https://github.com/lgragert/virtual-crossmatch.git)
5. Change directory to /conversion_tool/
6. Run command to start web server  
    - python3 manage.py runserver 8080  
7. Go to http://127.0.0.1:8080/ with your web browser. 

```






[GRAGERT LAB](https://hla.tulane.edu)