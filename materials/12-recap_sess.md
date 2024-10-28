---
title: Recap session (optional)
---

## Presentation and demo on basic command line usage

During the presentation, you will refresh your knowledge on command line usage. You will learn about Unix command line commands that you need to navigate in the filesystem, create and remove files, edit text documents. You will also learn about a `conda` like package and environment management solution (`micromamba`) that helps in software installations.

The presentation is accessible on Google Slides on the [following link](https://docs.google.com/presentation/d/1DTpToLWu4Q11-g2VwuhzqeWNEBYpZFgA0ze7fkSm_uE/edit?usp=sharing).

A short overview video from the Bioinformatics training Centre Team on Unix command line usage can be found on the [following link](https://drive.google.com/file/d/1qo_eNXaYSphU_Xpk0jzDHp1wJzyfSwAT/view)

::: {.callout-tip}
#### Most important command line commands

- `ls`: list folder contents
- `touch`: create an empty file (useful to check write permission)
- `nano`: open text editor
- `less`: view contents of a file (exit from the viewer by pressing `q`)
- `cp`: copy files and folders
- `mv`: move files and folders
- `grep`: search for pattern (e.g., word) in a text file
- `cd`: change directory
- `rm`: remove files and folders
- `cut`: extract sections of a file (e.g., a column of a table)
- `pwd`: print working directory (useful to see the full path)
:::

## Next-generation sequencing

During the metagenomics data analysis we perform several quality control steps and various types of trimmings and duplication removals. To understand why do you need to do these steps and understand the potential changes / biases you may expect, it is essential to understand the basics of Next-generation sequencing.

To refresh your knowledge on Next-generation sequencing techniques, please watch the tutorial video on the [following link](https://www.youtube.com/watch?v=mI0Fo9kaWqo)


## Summary

::: {.callout-tip}
#### Key Points

- Majority of the tools we use in metagenomics are command line tools, so a good working knowledge in Unix command line is obligatory.
- While is is more comfortable to use a text editor with graphical user interface (e.g., [Visual Studio Code](https://code.visualstudio.com)), the knowledge in pure command line based editor (e.g., `nano`) will come handy in many occasions.
- A clean and well-organised desktop or office makes the work more efficient, the same way a well-organised computational space makes our metagenomics analysis easier and more efficient. Package managers (like `conda`, `miniconda`, `mamba`, `micromamba`) combined with virtual environments help efficient software usage.
- Majority of the metagenomics data that we are analysing with bioinformatics tools come from Next-generation sequencing experiments. Knowing the fundamentals of these techniques significantly helps in the usage and understanding the potential benefits and caveats of several steps during the data analysis.
:::