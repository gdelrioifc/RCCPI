import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.sql.Time;
import java.util.ArrayList;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author Fonty and gdelrio
 */



public class GetRCC
{
 public ArrayList rcc;
 public GetRCC(String[] args)
  {
   try
    {
     Parameters.LoadParameters(args);

     if(!Parameters.loadCorrect)
      {
       return;
      }
     else
      {   
       if(Parameters.fileList != null)
        {
         int contador=0;
         BufferedReader fileList = new BufferedReader(new FileReader(Parameters.fileList));
         String line = "", clique_file="";
         while( (line = fileList.readLine())!= null)
          {
           String fields[] = line.split(",");
           Parameters.pdbFileName = fields[0];
           Parameters.chain = fields[1];
           Parameters.min=Double.parseDouble(fields[2]);
           Parameters.max=Double.parseDouble(fields[3]);
           if(fields[4].equalsIgnoreCase("yes")) Parameters.flagLateralChain=true;
           else if(fields[4].equalsIgnoreCase("no")) Parameters.flagLateralChain=false;
               
// La siguiente linea es la que uso Fernando para leer archivos compresos bajados del PDB
//                String localFilePath = Parameters.pathPdb + Parameters.pdbFileName.substring(1, 3).toLowerCase()+"/pdb"+Parameters.pdbFileName.toLowerCase()+".ent.gz";
           String localFilePath = Parameters.pathPdb + Parameters.pdbFileName;
           try
            {
             Peptide p = new Peptide( localFilePath, Parameters.chain, Parameters.flagLateralChain);
             if(p.residues.size() == 0)
              {
               System.out.println("file " + Parameters.pdbFileName + "  chain " + Parameters.chain + " is empty");
               continue;
              }
             contador=contador+1;
             Parameters.graphFileName = fields[0]+fields[1]+"-graph_"+String.valueOf(contador)+".txt";
             clique_file=fields[0]+fields[1]+"-cliques_"+String.valueOf(contador)+".txt";

             p.buildGraph(Parameters.min, Parameters.max, Parameters.distance_criterium, Parameters.pairs, Parameters.flagLateralChain, Parameters.graphFileName);
             Rcc r = new Rcc(Parameters.graphFileName,"tomita",clique_file);
//             System.out.println(line+","+r.printable);
            }
           catch (java.io.FileNotFoundException e)
            {
             System.out.println("Not found: " + fields[0] + "\tat\t"+localFilePath);
            }
           finally
            {
             continue;
            }
          }// end while
         fileList.close();
        }
       else
        {
         Peptide p = new Peptide (Parameters.pdbFileName, Parameters.chain, Parameters.flagLateralChain);

         if(p.residues.size() == 0)
          {
           System.out.println("file " + Parameters.pdbFileName + "  chain " + Parameters.chain + " is empty");
          }

         p.buildGraph(Parameters.min, Parameters.max, Parameters.distance_criterium, Parameters.pairs, Parameters.flagLateralChain, Parameters.graphFileName);
         Rcc r = new Rcc(Parameters.graphFileName,"tomita",Parameters.pdbFileName+"cliques.txt");
 
//         System.out.println(r.printable);
         rcc=r.getArrayList();
        
         return;
        }
      }
    }
   catch(Exception e)
    {
     System.err.println("Error at GetRCC constructor:");
     e.printStackTrace();
     return;
    }  
  }// end public GetRCC(String[] args) 
}
