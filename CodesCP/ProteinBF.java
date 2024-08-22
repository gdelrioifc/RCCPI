import java.io.*;
import java.util.*;

public class ProteinBF
{
 public int peptideStart, peptideEnd;
 private String pdbID, chain, pdbFileName;
 private ArrayList<ResidueBF> residues;

 public ProteinBF()
  {
   this.residues = new ArrayList<>();
  }

 public void setPDBFileName(String f)
  {
   this.pdbFileName=f;
  }

 public String getPDBFileName()
  {
   return this.pdbFileName;
  }

 public void setPDBID(String id)
  {
   this.pdbID=id;
  }

 public String getPDBID()
  {
   return this.pdbID;
  }

 public void setChain(String c)
  {
   this.chain=c;
  }

 public String getChain()
  {
   return this.chain;
  }

 public int getProteinSize()
  {
   return this.residues.size();
  }

 public void addResidue(ResidueBF r)
  {
   this.residues.add(r);
  }

 public ResidueBF getResidue(int i)
  {
   return (ResidueBF)this.residues.get(i);
  }

 public ArrayList<ResidueBF> getResidues()
  {
   return this.residues;
  }

 public String[] getPeptideFile(int pep_ini, int pep_len)
  {
   String[] result=new String[3];
   try
    {
     int i=0, count=0;
     String line="", resN="", resS="";
     ProteinBF peptide = this.getPeptide(pep_ini, pep_len);
     ResidueBF r1=null, r2=null;
     ArrayList<ResidueBF> p_residues=peptide.getResidues();
     BufferedReader infile = new BufferedReader(new FileReader(this.pdbFileName));
     BufferedWriter outfile = null;
     result[0]="peptide_"+String.valueOf(count)+".pdb";
     File outn = null;
     outn = new File(result[0]);
     while(outn.exists())
      {
       count=count+1;
       result[0]="peptide_"+String.valueOf(count)+".pdb";
       outn = new File(result[0]);
      }
     result[1]=((ResidueBF)p_residues.get(0)).toString();
     result[2]=((ResidueBF)p_residues.get(p_residues.size()-1)).toString();
     outfile= new BufferedWriter(new FileWriter(result[0]));

     while((line=infile.readLine())!=null)
      {
       if(line.startsWith("ATOM"))
        {
         chain=line.substring(21,22).trim();
         if(chain.equals(this.chain))
          {
           resN=line.substring(17,20).trim();
           resS=line.substring(22,26).trim();
           r1=new ResidueBF(resN,resS,chain);
           breakable:
           for(i=0; i<p_residues.size();i++)
            {
             r2=(ResidueBF)p_residues.get(i);
             if(r1.compareResidueNames(r2))
              {
               outfile.write(line+"\n");
               break breakable;
              }
            }
          }
        }
      }
     infile.close();
     outfile.flush();
     outfile.close();
    }
   catch(Exception e)
    {
     System.err.println("Error @ ProteinBF.getPeptideFile:");
     e.printStackTrace();
    }
   return result;
  }// end getPeptideFile

 public ProteinBF getPeptide(int pep_ini, int pep_len)
  {
   ProteinBF result = new ProteinBF();
   result.setPDBFileName(this.pdbFileName);
   result.setPDBID(this.pdbID);
   result.setChain(this.chain);
   ResidueBF r=null;
   try
    {
     String atom="", coord="";
     if(this.residues!=null)
      {
//System.err.println("Reading protein to get peptide:"+this.toString());
       int i=0;
       result.peptideStart=pep_ini;
       if((pep_ini+pep_len)>this.residues.size()) result.peptideEnd=this.residues.size();
       else result.peptideEnd=pep_ini+pep_len;

       for(i=pep_ini; i<result.peptideEnd; i++)
        {
         r=(ResidueBF)this.residues.get(i);
         result.addResidue(r);
//System.err.println(r.toString());
        }
//System.err.println("Generated peptide:"+result.toString()+". Peptide's atoms:"+result.getProteinSize());
      }
    }
   catch(Exception e)
    {
     System.err.println("Error @ ProteinBF.getPeptide:");
     e.printStackTrace();
    }
   return result;
  }

 public ProteinBF setPeptide(int pep_ini, int pep_len)
  {
   ProteinBF result = new ProteinBF();

   result.setPDBFileName(this.pdbFileName);
   result.setPDBID(this.pdbID);
   result.setChain(this.chain);
   result.peptideStart=pep_ini;
   if((pep_ini+pep_len)>this.residues.size()) result.peptideEnd=this.residues.size();
   else result.peptideEnd=pep_ini+pep_len;

   return result;
  }

 public int getPeptideStart()
  {
   return this.peptideStart;
  }

 public int getPeptideEnd()
  {
   return this.peptideEnd;
  }

 public ProteinBF duplicate()
  {
   ProteinBF result = this;

   return result;
  }

 public String getCADistance(ProteinBF pep)
  {
   String result="";
   try
    {
     int i=0, j=0;
     double x1=0.0, y1=0.0, z1=0.0, x2=0.0, y2=0.0, z2=0.0, distance=0.0, smallest=1000000.0;
     String coord1="", coord2="";
     String[] data1=null, data2=null;
     ResidueBF r1=null, r2=null;

     for(i=0; i<this.residues.size(); i++)
      {
       r1=this.getResidue(i);
       coord1=(String)r1.getCoordsCA(r1);
       data1=coord1.split(",");
       x1=Double.parseDouble(data1[0]);
       y1=Double.parseDouble(data1[1]);
       z1=Double.parseDouble(data1[2]);


       for(j=0; j<pep.residues.size(); j++)
        {
         r2=pep.getResidue(j);
         coord2=(String)r2.getCoordsCA(r2);
         data2=coord2.split(",");
         x2=Double.parseDouble(data2[0]);
         y2=Double.parseDouble(data2[1]);
         z2=Double.parseDouble(data2[2]);


         distance=Math.pow((x1-x2),2.0)+Math.pow((y1-y2),2.0)+Math.pow((z1-z2),2.0);
         distance=Math.sqrt(distance);

         if(distance<smallest)
          {
           smallest=distance;
           result=r1.getName()+r1.getResSeq()+r1.getChain()+"-"+r2.getName()+r2.getResSeq()+r2.getChain()+","+smallest;
          }
        }
      }
    }
   catch(Exception e)
    {
     System.err.println("Error @ ProteinBF.getCADistance:");
     e.printStackTrace();
    }
   return result;
  }// end getCADistance

 public String getSmallestDistance(ProteinBF pep)
  {
   String result="";
   try
    {
     int i=0, j=0, k=0, l=0;
     double x1=0.0, y1=0.0, z1=0.0, x2=0.0, y2=0.0, z2=0.0, distance=0.0, smallest=1000000.0;
     String coord1="", coord2="";
     ArrayList coords1=null, coords2=null;
     String[] data1=null, data2=null;
     ResidueBF r1=null, r2=null;

     for(i=0; i<this.residues.size(); i++)
      {
       r1=this.getResidue(i);
       coords1=r1.getCoords();
       for(j=0; j<coords1.size(); j++)
        {
         data1=((String)coords1.get(j)).split(",");
         x1=Double.parseDouble(data1[0]);
         y1=Double.parseDouble(data1[1]);
         z1=Double.parseDouble(data1[2]);

         for(k=0; k<pep.residues.size(); k++)
          {
           r2=pep.getResidue(k);
           coords2=r2.getCoords();
           for(l=0; l<coords2.size(); l++)
            {
             data2=((String)coords2.get(l)).split(",");
             x2=Double.parseDouble(data2[0]);
             y2=Double.parseDouble(data2[1]);
             z2=Double.parseDouble(data2[2]);


             distance=Math.pow((x1-x2),2.0)+Math.pow((y1-y2),2.0)+Math.pow((z1-z2),2.0);
             distance=Math.sqrt(distance);

             if(distance<smallest)
              {
               smallest=distance;
               result=r1.getName()+r1.getResSeq()+r1.getChain()+"-"+r2.getName()+r2.getResSeq()+r2.getChain()+","+smallest;
              }
            }
          }
        }
      }
    }
   catch(Exception e)
    {
     System.err.println("Error @ ProteinBF.getCADistance:");
     e.printStackTrace();
    }
   return result;
  }// end getSmallestDistance 

 public String getSmallestAvgDistance(ProteinBF pep)
  {
   String result="";
   try
    {
     int i=0, j=0, k=0, l=0;
     double x1=0.0, y1=0.0, z1=0.0, x2=0.0, y2=0.0, z2=0.0, distance=0.0, avg=0.0, smallest=1000000000.0;
     String coord1="", coord2="";
     ArrayList coords1=null, coords2=null;
     String[] data1=null, data2=null;
     ResidueBF r1=null, r2=null;

     for(i=0; i<this.residues.size(); i++)
      {
       r1=this.getResidue(i);
       coords1=r1.getCoords();
       for(j=0; j<coords1.size(); j++)
        {
         data1=((String)coords1.get(j)).split(",");
         x1=Double.parseDouble(data1[0]);
         y1=Double.parseDouble(data1[1]);
         z1=Double.parseDouble(data1[2]);

         avg=0;
         for(k=0; k<pep.residues.size(); k++)
          {
           r2=pep.getResidue(k);
           coords2=r2.getCoords();
           for(l=0; l<coords2.size(); l++)
            {
             data2=((String)coords2.get(l)).split(",");
             x2=Double.parseDouble(data2[0]);
             y2=Double.parseDouble(data2[1]);
             z2=Double.parseDouble(data2[2]);

             distance=Math.pow((x1-x2),2.0)+Math.pow((y1-y2),2.0)+Math.pow((z1-z2),2.0);
             distance=Math.sqrt(distance);

             avg=avg+distance;
            }
          }
         avg=avg/(1.0*pep.residues.size());
         if(avg<smallest)
          {
           smallest=avg;
           result=String.valueOf(smallest)+","+((ResidueBF)pep.residues.get(0)).toString()+"-"+((ResidueBF)pep.residues.get(pep.residues.size()-1)).toString();
          }
        }
      }
    }
   catch(Exception e)
    {
     System.err.println("Error @ ProteinBF.getSmallestAvgDistance:");
     e.printStackTrace();
    }
   return result;
  }// end getSmallestAvgDistance

 public String toString()
  {
   String result="";

   result=this.pdbID+","+this.chain+","+((ResidueBF)this.residues.get(0)).getName()+((ResidueBF)this.residues.get(0)).getResSeq()+((ResidueBF)this.residues.get(0)).getChain()+","+((ResidueBF)this.residues.get(this.residues.size()-1)).getName()+((ResidueBF)this.residues.get(this.residues.size()-1)).getResSeq()+((ResidueBF)this.residues.get(this.residues.size()-1)).getChain();

   return result;
  }
}
