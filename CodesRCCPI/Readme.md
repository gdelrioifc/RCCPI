# CodesRCCPI
<h2>This directory contains the Java files required to identify peptides that match protein-protein interfaces.</h2>

<p>There are two types of files in this dir: java and jar files. The first ones you'll need them to compile an executable Java file; the second ones you'll need to compile and execute the Java file as explained below.</p>

<p>After downloading all these files in your computer, go to the location where you save the files and type in your terminal:</p>

<p>javac -cp weka.jar:autoweka.jar:. PrPrInterfacePredictorBFRCCFromFile.java</p>

<p>Depending on your Java version, you may get some warning messages; you may ignore them.</p>

<p>That command line will generate a .class file for each of your .java files.</p>

<p>To execute the code useful to predict the peptides from one protein that interact with another protein, you may type in your terminal:</p>

<p>java -cp weka.jar:autoweka.jar:. PrPrInterfacePredictorBFRCCFromFile</p>

<p>That will provide you with the instructions to compute those peptides.</p>
