import java.io.*;
import java.util.*;

import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.*;

import weka.classifiers.meta.AutoWEKAClassifier;
import weka.core.*;
import weka.core.converters.ConverterUtils.DataSource;
import weka.core.converters.ConverterUtils.DataSink;

public class PrPrInterfacePredictorBFRCCFromFile
{
 public static void main(String[] args) throws Exception
  {
   if(args.length==1)
    {
     int i=0, processors = Runtime.getRuntime().availableProcessors();
     long waiting_time=2;
     String line="", result="";
     String[] token=null;
     GetMaximumScoreRCCRunnable gmsrcc=null;
     ExecutorService executor=Executors.newFixedThreadPool(processors);;
     BufferedReader infile = new BufferedReader(new FileReader(args[0]));

     while((line=infile.readLine())!=null)
      {
       if(line.indexOf(",")>-1)
        {
         token=line.split(",");
         if(token.length==11)
          {
           i=i+1;
           if(i>0 && i<=processors)
            {
             gmsrcc=new GetMaximumScoreRCCRunnable(token);
             executor.execute(gmsrcc);
             if(i==processors)
              {
               executor.awaitTermination(waiting_time,TimeUnit.MINUTES);
               i=0;
              }
            }
          }
        }
      }
     infile.close();
     executor.awaitTermination(waiting_time,TimeUnit.MINUTES);
//     executor.shutdown();
    }
   else
    {
     System.err.println("Usage:\njava PrPrInterfacePredictorBFRCCFromFile <infile>");
     System.err.println("<infile> a CSV file containing the following info in the specified order: <pdb>,<chain1>,<chain2>,<pep_len>,<hits>,<weka_model>,<lateralSideChain>,<minDist>,<maxDist>,<operation>,<outfile>.");
     System.err.println("<pdb> the name of a file in PDB format; e.g., 5dnk.pdb (Please include the full path where this file is located in your computer).");
     System.err.println("<chain1> a letter specifying the chain in <pdb> part of a protein-protein interaction (PPI); e.g., A.");
     System.err.println("<chain2> a letter specifying the other chain in <pdb> part of a PPI; e.g., B.");
     System.err.println("<pep_len> intereger specifying the peptide length to be generated; e.g., 20.");
     System.err.println("<hits> integer specifying how many peptides clser to <cutoff> to report; e.g., 3.");
     System.err.println("<weka_model> file containing a model trained with AutoWeka; e.g.");
     System.err.println("<lateralSideChain> yes or no. This parameter should correspond with the specified <weka_model>.");
     System.err.println("<minDist> interger, usually 0. This parameter should correspond with the specified <weka_model>.");
     System.err.println("<maxDist> integer. This parameter should correspond with the specified <weka_model>.");
     System.err.println("<operation> sum or concat, according to the specified <weka_model>.");
     System.err.println("<outfile> name of the file where to save the results.");
     System.err.println("For each line, this program will launch a process and will do that for any subsequent line, untile all the processors");
     System.err.println("In the local machine are busy. Then, it will wait for the processes to finish, and repeat the same process until all");
     System.err.println("all the lines are executed.");
     System.err.println("The program implements a brute force algorithm to identify the region (peptides) that have the largest");
     System.err.println("score for peptide-protein prediction based on RCC; such region may be the interface between <chain1>-<chain2>.");
     System.err.println("The output will include the peptides from each chain in the pdb that is predicted");
     System.err.println("to include the PPI interface. The length of the peptides will range from <pep_min> and <pep_max>.");
     System.err.println("");
     System.err.println("");
    }
  }
}
