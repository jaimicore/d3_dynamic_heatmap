################
## D3 heatmap ##
################
##
## This script produces an interactive heatmap in D3 from a numeric table.
##
## Input:
##	Mandatory:
##		- A tab-delimited file (input.tab.file)
##	Optional:
##		- Title displayed at the top of the page (page.title)
##		- The d3 base path, this could be a file or a url (d3.base.file)
##		- The results directory (results.dir)
##		- The path to the HTML template (heatmap.template.file)
##		- The color palette (color.palette). The same for the colorBrewer palettes (http://colorbrewer2.org). Default: YlOrRd
##		- The number of color classes (color.palette). The same for the colorBrewer palettes (http://colorbrewer2.org). Default: 9
##
## How to run it (assuming you run it in the repository directory):
cat R/Matrix_to_dynamic_d3_heatmap.R | 
/usr/bin/R --slave --no-save --no-restore --no-environ --args " \
results.dir = 'D3_demo_results' ; \
input.tab.file = 'Demo_files/RNAseq_example.tab' ; "
##
## Output: a html file with the dynamic D3 heatmap. Use firefox to visualize it.


