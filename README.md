# MarineOmics Website
Here you will find a growing and evolving collection of guidelines, tutorials, recommended readings, and method evaluations for various ’omics analyses that are commonly used in studies of nonmodel systems (particularly marine species). Instead of providing rigid “best practices”, we aim to encourage data exploration and rigorous quality control through “best principles”. The guidelines are written by members of the MarineOmics working group or the community, then shared for peer review. Comments on current website content or suggestions for new content can be made through the [Discussion Forum](https://github.com/MarineOmics/marineomics.github.io/discussions).

## Contributing to the website
The website design and editing workflow is inspired by GitHub sites managed by the [Grunwald Lab](https://github.com/grunwaldlab/Population_Genetics_in_R). To add a new page of content, an RMarkdown file is created and then Knitted into an HTML file. This allows for the integration of code (R, Python, Bash) with Markdown formatting. The `_site.yml` file controls the style of the website and the Navigation Bar structure. To add your new page to the website, you must add the HTML file name into the `_site.yml` file.  

If your are a member of the MarineOmics GitHub organization, you can contribute to the website repository as usual:  
```
git clone https://github.com/MarineOmics/marineomics.github.io.git
# Make changes 
git add [.Rmf file and .html file]
git commit ...
git push ...
```

### Including images  
To include images into your page, add the image file to the `images/` folder and reference the image in your RMarkdown file.
```
![Figure comparing the distribution of sequencing reads mapped to a reference genome.](images/mec16077-fig-0001-m.jpg)
```
### Including code
If your code requires an input file, you should create a new folder in the GitHub repo to hold that data file. Try to do this sparingly, and be aware of file size restrictions.  

When you Knit an RMarkdown file and it generates figures, the plots for these figures should be output into a new folder with `_files` in the name. Make sure to add this folder to your Github commit.  

### Including references  
You can add Bibtex references to the file `common-bib_01.bib`, and then refer to them in your RMarkdown by adding `bibliography: common-bib_01` to the header and citing them in the RMarkdown document. Check out existing pages for examples.  
