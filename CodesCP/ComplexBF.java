import java.io.*;
import java.util.*;

public class ComplexBF
{
 private String pdbID, pdbFileName;

 public ComplexBF(String f)
  {
   if(f.indexOf(File.separator)>-1) pdbID=f.substring(f.indexOf(File.separator)+1,f.lastIndexOf("."));
   else pdbID=f.substring(0,f.lastIndexOf("."));
   pdbFileName=f;
  }

 public ProteinBF setProteinBF(String c)
  {
   ProteinBF result = new ProteinBF();

   result.setPDBFileName(this.pdbFileName);
   result.setPDBID(this.pdbID);
   result.setChain(c);

   return result;
  }// end setProteinBF

 public ProteinBF getProteinBF(String c)
  {
   ProteinBF result = new ProteinBF();
   result.setPDBID(this.pdbID);
   result.setChain(c);
   try
    {
     String line="", chain="", resN="", resS="";
     ResidueBF r=null;
     BufferedReader infile = new BufferedReader(new FileReader(pdbFileName));

     while((line=infile.readLine())!=null)
      {
       if(line.startsWith("ATOM"))
        {
         chain=line.substring(21,22).trim();
         if(chain.equals(c))
          {
           resN=line.substring(17,20).trim();
           resS=line.substring(22,26).trim();
           if(r==null) 
            {
             r=new ResidueBF(resN,resS,chain);
             r.add(line);
             result.addResidue(r);
            }
           else
            {
             if(!r.add(line))
              {
               r=new ResidueBF(resN,resS,chain);
               r.add(line);
               result.addResidue(r);
              }
            }
          }
        }
      }
     infile.close();
     result.peptideStart=0;
     result.peptideEnd=result.getProteinSize()-1;
    }
   catch(Exception e)
    {
     System.err.println("Error @ Complex.getProteinBF:");
     e.printStackTrace();
    }
   return result;   
  }
}
