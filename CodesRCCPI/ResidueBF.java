import java.io.*;
import java.util.*;

public class ResidueBF
{
 private String resName="", resSeq="", chain="";
 public ArrayList atoms, coords;

 public ResidueBF(String resN, String resS, String c)
  {
   atoms=new ArrayList();
   coords=new ArrayList();
   this.resName=resN;
   this.resSeq=resS;
   this.chain=c;
  }

 public String getName()
  {
   return this.resName;
  }

 public String getResSeq()
  {
   return this.resSeq;
  }

 public String getChain()
  {
   return this.chain;
  }

 public String toString()
  {
   String result=this.resName+this.resSeq+"_"+this.chain;
   return result;
  }

 public boolean compareResidueNames(ResidueBF r)
  {
   boolean result=false;
   if(this.getName().equals(r.getName()) && this.getResSeq().equals(r.getResSeq()) && this.getChain().equals(r.getChain())) result=true;
   return result;
  }

 public boolean add(String line)
  {
   boolean result = false;
   try
    {
     String resN=line.substring(17,20).trim(), resS=line.substring(22,26).trim();
     if(this.resSeq.equals(resS) && this.resName.equals(resN))
      {
       result=true;

       String atom=line.substring(12,16).trim();
       String coord=line.substring(30,38).trim()+","+line.substring(38,46).trim()+","+line.substring(46,54).trim();
       atoms.add(atom);
       coords.add(coord);
      }
    }
   catch(Exception e)
    {
     System.err.println("Error @ ResidueBF.add:");
     e.printStackTrace();
    }
   return result;
  }

 public String getCoordsCA(ResidueBF r)
  {
   String result="";

   int i=0;
   for(i=0; i<r.atoms.size(); i++)
    {
     if(((String)r.atoms.get(i)).equals("CA"))
      {
       result=(String)r.coords.get(i);
      }
    }

   return result;
  }

 public void setCoords(ArrayList c)
  {
   this.coords=c;
  }

 public ArrayList getCoords()
  {
   return this.coords;
  }// end getCoords

 public void setAtoms(ArrayList a)
  {
   this.atoms=a;
  }

 public int countAtoms()
  {
   return this.atoms.size();
  }
}
