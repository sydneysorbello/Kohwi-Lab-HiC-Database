# Kohwi-Lab-HiC-Database
This repository provides code for the Kohwi Lab HiC Database. This database was produced for an assignment in Boston University's graduate level course BF768. The database was created to support the Kohwi Lab of Columbia University. No experimental data is provided, but the flask app and html that supported the host website are available here.

The database was produced by 4 students in Boston University's Master in Bioinformatics Program: Sydney Sorbello, Katie Kitrick, Dhruvi Joshi, and Sofiya Patra. 

## Application Tabs / Features

The web application consists of five main tabs:

Welcome — Overview of the database purpose, lab information, and developer team.

Hi-C — Direct interaction with the SQL database storing Hi-C interaction data. Users can submit queries, download resulting data as CSV, visualize contact matrices, and generate interaction frequency plots.

Genome — Embedded UCSC Genome Browser for Drosophila genome visualization and exploration.

Help — Instructions on how to use the database, run queries, and interpret outputs.

Login / Admin — Authenticated upload/edit interface for lab members to add or modify experimental datasets.

## Home Page Example

Below is a screenshot of the website Home Page:

![Home Page Screenshot](https://github.com/sydneysorbello/Kohwi-Lab-HiC-Database/blob/main/app/static/plots/Home_Page_Screenshot.png?raw=true)

## SQL Query Interface

Database query is made easy through drop down boxes containing sample and genomic location selection:

![SQL_Query_Interface](https://github.com/sydneysorbello/Kohwi-Lab-HiC-Database/blob/main/app/static/plots/SQL_Query_Interface.png?raw=true)

## Visualization of Hi-C Matrix

The Generate Plot feature seemlessly visualizes a Hi-C Matrix given the input parameters:

![HiC_Matrix_Example](https://github.com/sydneysorbello/Kohwi-Lab-HiC-Database/blob/main/app/static/plots/Generate_Plot_Example.png?raw=true)

## Video Tutorial

Walkthrough video tutorial available here:

https://github.com/sydneysorbello/Kohwi-Lab-HiC-Database/tree/main/website_screenshots
