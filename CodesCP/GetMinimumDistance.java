import java.io.*;
import java.util.*;


class GetMinimumDistance implements Runnable
{
 String[] args;

 public GetMinimumDistance(String[] a)
  {
   args=a;
  }

 public static void lowestDistance(ProteinBF p1, ProteinBF p2, int pep_len, int hits, String outName)
  {
   try
    {
     boolean valid=true;
     int p1_size=p1.getProteinSize(), p2_size=p2.getProteinSize(), pep_start=0, pep_end=0, pep_ini=0, i=0;
     double d1=0.0, distance=10000.0, smallest=distance;
     Double score=null;
     String[] data1=null;
     List<Double> tmp=null;
     Iterator<Double> it = null;
     Hashtable result = null;
     ProteinBF pep_current=null;
     BufferedWriter outfile=new BufferedWriter(new FileWriter(outName));

// obtaining best interacting peptide from protein 1 to protein 2
     result=new Hashtable();
     pep_ini=0;
     while(valid)
      {
       // iterate over all possible peptides of protein 1
       pep_current=p1.getPeptide(pep_ini, pep_len);

       if(p1.getProteinSize() < (pep_ini+pep_len)) valid=false;
       else
        {
         data1=p2.getSmallestAvgDistance(pep_current).split(",");
         d1=Double.parseDouble(data1[0]);
         if(d1<smallest) result.put(new Double(d1),data1[1]);
        }
       pep_ini=pep_ini+1;
      }

     tmp = Collections.list(result.keys());
     Collections.sort(tmp);
     it = tmp.iterator();
     i=0;
     breakable:
     while(it.hasNext())
      {
       i=i+1;
       if(i<=hits)
        {
         score=(Double)it.next();
         outfile.write(score.toString()+":"+(String)result.get(score)+"\n");
        }
       else break breakable;
      }
     outfile.flush();
     outfile.close();
/*
// obtaining best interacting peptide from protein 2 to protein 1
     valid=true;
     pep_ini=0;
     distance=10000.0;
     smallest=distance;
     pep_current=null;
     result=new Hashtable();
     while(valid)
      {
       pep_current=p2.getPeptide(pep_ini, pep_len);
       if(p1.getProteinSize() < (pep_ini+pep_len)) valid=false;
       else
        {
         data1=p1.getSmallestAvgDistance(pep_current).split(",");
         d1=Double.parseDouble(data1[0]);
         if(d1<smallest) result.put(new Double(data1[0]),data1[1]);
        }
       pep_ini=pep_ini+1;
      }

     tmp = Collections.list(result.keys());
     Collections.sort(tmp);
     it = tmp.iterator();
     i=0;
     breakable:
     while(it.hasNext())
      {
       i=i+1;
       if(i<=hits)
        {
         score=(Double)it.next();
         System.out.println(score.toString()+":"+(String)result.get(score));
        }
       else break breakable;
      }
*/
    }
   catch(Exception e)
    {
     System.err.println("Error @ PrPrInterfacePredictor.lowestDistance:");
     e.printStackTrace();
    }
  }// end lowestDistance
 
 public void run()
  {
   try
    {
     if(args.length==6)
      {
       ComplexBF complex = new ComplexBF(args[0]);
       int pep_len=Integer.parseInt(args[3]), hits=Integer.parseInt(args[4]);

       ProteinBF p1 = complex.getProteinBF(args[1]), p2 = complex.getProteinBF(args[2]);
//System.err.println("P1 data:"+p1.toString());
//System.err.println("P2 data:"+p2.toString());
       lowestDistance(p1,p2,pep_len,hits, args[5]);
      }
    }
   catch(Exception e)
    {
     System.err.println("Error @ GetMinimumDistance.run:");
     e.printStackTrace();
    }
  }
}
