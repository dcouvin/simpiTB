# simpiTB
A pipeline for the analysis of whole genome sequences belonging to Mycobacterium tuberculosis complex

## Singularity
Use the following command to build the Singularity image from this current repository:
`sudo singularity build -F simpiTB.simg Singularity/simpiTB.def`

or 

`sudo singularity config fakeroot --add <user>`

`singularity build -f simpiTB.simg Singularity/simpiTB.def`
