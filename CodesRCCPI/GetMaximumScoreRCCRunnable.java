import java.io.*;
import java.util.*;

import weka.classifiers.meta.AutoWEKAClassifier;
import weka.core.*;
import weka.core.converters.ConverterUtils.DataSource;

class GetMaximumScoreRCCRunnable implements Runnable
{
 String[] args;

 public GetMaximumScoreRCCRunnable(String[] a)
  {
   args=a;
  }

 public static ArrayList getResultingRCC(ArrayList rcc1, ArrayList rcc2, String operation)
  {
   ArrayList result = new ArrayList();
   try
    {
     int i=0;
     int d1=0, d2=0, d=0;

     if(operation.equalsIgnoreCase("concat"))
      {
       for(i=0; i<rcc1.size(); i++) result.add(String.valueOf(((Integer)rcc1.get(i)).intValue()));
       for(i=0; i<rcc2.size(); i++) result.add(String.valueOf(((Integer)rcc2.get(i)).intValue()));
      }
     else if(operation.equalsIgnoreCase("sum"))
      {
       for(i=0; i<rcc1.size(); i++) 
        {
         d1=((Integer)rcc1.get(i)).intValue();
         d2=((Integer)rcc2.get(i)).intValue();
         d=d1+d2;
         result.add(String.valueOf(d));
        }
      }
    }
   catch(Exception e)
    {
     System.err.println("Error @ GetMaximumScoreRCCRunnable.getResultingRCC:");
     e.printStackTrace();
    }
   return result;
  }// end getResultingRCC

 public static String createARFFFile(ArrayList rcc)
  {
   String result="";
   try
    {
     String outNameTag="testFile_", outName="";
     File file = null;
     BufferedWriter outfile=null;
     ArrayList<Attribute> atts = new ArrayList<Attribute>();
     ArrayList<String>    attVals = new ArrayList<String>();
     Instances            data;
     double[]             vals;
     int                  i, count=0;

     outName=outNameTag+String.valueOf(count)+".arff";
     file=new File(outName);
     while(file.exists())
      {
       count=count+1;
       outName=outNameTag+String.valueOf(count)+".arff";
       file=new File(outName);
      }

     if(rcc.size()==26)
      {
       for(i=1;i<27;i++) atts.add(new Attribute("RRC"+String.valueOf(i))); 
      }
     else if(rcc.size()==52)
      {
       for(i=1;i<53;i++) atts.add(new Attribute("RRC"+String.valueOf(i)));
      }
     attVals.add("P"); attVals.add("N");
     atts.add(new Attribute("Class", attVals));
     
     data = new Instances("ProtPepRelation", atts, 0);

     vals = new double[data.numAttributes()];

     if(rcc.size()==26)
      {
       for(i=0;i<26;i++) vals[i]=Double.parseDouble((String)rcc.get(i));
      }
     else if(rcc.size()==52)
      {
       for(i=0;i<52;i++) vals[i]=Double.parseDouble((String)rcc.get(i));
      }
     data.add(new DenseInstance(1.0, vals));

     outfile = new BufferedWriter(new FileWriter(outName));
     outfile.write(data.toString());
     outfile.flush();
     outfile.close();
     result=outName;
    }
   catch(Exception e)
    {
     System.err.println("Error @ GetMaximumScoreRCCRunnable.createARFFFile:");
     e.printStackTrace();
    }
   return result;
  }// end createARFFFile

 public static double getPrediction(String fModel, String arff)
  {
   double result=-1.0;
   try
    {
     int i=0, j=0;
     double value=0.0;
     double[] predictVal=null;
     String distribution="", prediction="";
     Instance instance=null;
     AutoWEKAClassifier model = (AutoWEKAClassifier)weka.core.SerializationHelper.read(fModel);
     DataSource source= new DataSource(arff);
     Instances Testset=source.getDataSet();
     Testset.setClassIndex(Testset.numAttributes()-1);

     for(i=0; i<Testset.numInstances(); i++)
      {
       instance=Testset.instance(i);
       value=model.classifyInstance(instance);
       predictVal=model.distributionForInstance(instance);
       prediction=Testset.classAttribute().value((int)value);
       result=predictVal[0];
      }
    }
   catch(Exception e)
    {
     System.err.println("Error @ GetMaximumScoreRCCRunnable.getPrediction:");
     e.printStackTrace();
    }
   return result;
  }// end getPrediction

 public static boolean isNotOnly0s(ArrayList rcc)
  {
   boolean result=true;
   try
    {
     int i=0, count=0;
     String v="";

     for(i=0; i<rcc.size();i++)
      {
       v=String.valueOf(((Integer)rcc.get(i)).intValue()).trim();
       if(v.equals("0")) count=count+1;
      }
     if(count==26) result=false;
    }
   catch(Exception e)
    {
     System.err.println("Error @ PrPrInterfacePredictorBF.isNotOnly0s:");
     e.printStackTrace();
    }
   return result;
  }// end isNotOnly0s

 public static Hashtable maximumScoreRCC(String pdb, String chain1, String chain2, int pep_len, String model, String lateral, String minDist, String maxDist, String operation)
  {
   Hashtable result = new Hashtable();
   try
    {
     boolean valid=true;
     int pep_ini=0;
     double score=0.0, largest=-1.0;
     String arffName="";
     String[] args=null, pep_currentData=null;
     File deleteFile=null;
     GetRCC getRCC=null;
     ArrayList rcc=null, rcc1=null, rcc2=null, peptides=null;
     ComplexBF complex = new ComplexBF(pdb);
     ProteinBF p1 = complex.getProteinBF(chain1), p2 = complex.getProteinBF(chain2);
     p1.setPDBFileName(pdb); 
     p1.setChain(chain1);
     p2.setPDBFileName(pdb);
     p2.setChain(chain2);

// obtaining best interacting peptide from protein 1 to protein 2
     args=new String[]{"-pdb",pdb,"-chain",chain2,"-pathTomita","/mnt/sde/gdelrio/quick-cliques/bin/","-lateralChain",lateral,"-minDist",minDist,"-maxDist",maxDist};
     Parameters.LoadParameters(args);

     if(Parameters.loadCorrect)
      {
       getRCC = new GetRCC(args);
       rcc2=getRCC.rcc;
      }

     pep_ini=0;
     while(valid)
      {
       // iterate over all possible peptides of protein 1
       if(p1.getProteinSize()>(pep_ini+pep_len))
        {
         pep_currentData=p1.getPeptideFile(pep_ini, pep_len);
         deleteFile=new File(pep_currentData[0]);
         args=new String[]{"-pdb",pep_currentData[0],"-chain",chain1,"-pathTomita","/mnt/sde/gdelrio/quick-cliques/bin/","-lateralChain",lateral,"-minDist",minDist,"-maxDist",maxDist};
         Parameters.LoadParameters(args);
         if(Parameters.loadCorrect)
          {
           getRCC = new GetRCC(args);
           rcc1=getRCC.rcc;
           deleteFile.delete();
           if(isNotOnly0s(rcc1))
            {
             rcc=getResultingRCC(rcc1,rcc2,operation);
             arffName=createARFFFile(rcc);
             score=getPrediction(model,arffName);
             if(score>=largest)  
              {
               largest=score;
               if(result.containsKey(new Double(largest)))
                {
                 peptides=(ArrayList)result.get(new Double(largest));
                 peptides.add(pep_currentData[1]+"-"+pep_currentData[2]);
                 result.put(new Double(largest),peptides);
                }
               else
                {
                 peptides=new ArrayList();
                 peptides.add(pep_currentData[1]+"-"+pep_currentData[2]);
                 result.put(new Double(largest),peptides);
                }
              }
             deleteFile=new File(arffName);
             deleteFile.delete();
            }
          }
         pep_ini=pep_ini+1;
        }
       else valid=false;
      }
    }
   catch(Exception e)
    {
     System.err.println("Error @ GetMaximumScoreRCCRunnable.maximumScoreRCC:");
     e.printStackTrace();
    }
   return result;
  }// end maximumScoreRCC

 public void run()
  {
   try
    {
     int i=0, j=0;
     Double score=null;
     String minDist=args[7], maxDist=args[8], operation=args[9], outfileName=args[10];
     String model=args[5], lateral=args[6], pdb=args[0], chain1=args[1], chain2=args[2];
     int pep_len=Integer.parseInt(args[3]), hits=Integer.parseInt(args[4]);
     Hashtable best_peptides=null;
     ArrayList peptides=null;
     List<Double> tmp=null;
     Iterator<Double> it=null;
     BufferedWriter outfile = new BufferedWriter(new FileWriter(outfileName));

     best_peptides=maximumScoreRCC(pdb, chain1, chain2, pep_len, model, lateral, minDist, maxDist, operation);

     if(best_peptides.size()>0)
      {
       tmp = Collections.list(best_peptides.keys());
       Collections.sort(tmp,Collections.reverseOrder());
       it = tmp.iterator();
       outfile.write("Different scores reported in "+args[0]+":"+best_peptides.size()+"\n");
       breakable:
       while(it.hasNext())
        {
         score=(Double)it.next();
         peptides=(ArrayList)best_peptides.get(score);
         for(i=(peptides.size()-1); i>=0; i--)
          {
           j=j+1;
           if(j<=hits)
            {
             outfile.write(score.toString()+":"+(String)peptides.get(i)+"\n");
            }
           else break breakable;
          }
        }
       outfile.flush();
       outfile.close();
      }
    }
   catch(Exception e)
    {
     System.err.println("Error @ GetMaximumScoreRCCRunnable.run:");
     e.printStackTrace();
    }
  }
}
