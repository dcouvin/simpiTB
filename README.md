# simpiTB
A pipeline for the analysis of whole genome sequences belonging to Mycobacterium tuberculosis complex
![simipTB](logo_simpiTB.png)

## Galaxy
A Galaxy wrapper has been developed to make the simpiTB tool more accessible. The tool is available at [Galaxy KaruBioNet](http://calamar.univ-ag.fr/c3i/galaxy_karubionet.html)

## Singularity
Use the following command to build the Singularity image from this current repository:
`sudo singularity build -F simpiTB.simg Singularity/simpiTB.def`

or 

`sudo singularity config fakeroot --add <user>`

`singularity build -f simpiTB.simg Singularity/simpiTB.def`
