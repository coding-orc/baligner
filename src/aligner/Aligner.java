/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package aligner;

import core.sb.SB;
import core.sb.Seq;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.lang.ProcessBuilder.Redirect;
import java.net.URL;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

/**
 *
 * @author John
 */
public class Aligner {

    /**
     * @param args the command line arguments
     */
    
    //with profile
    //-in /Users/John/Dropbox/code/Java/Aligner/input.txt -minp 6
    public static void main(String[] args) {
        try{
            
            header();
            String IN_FILE = "";
            int MIN_PROF_SIZE = 0;
            boolean DOPROFILE = false;

            for(int x=0; x<args.length; x++){
                if(args[x].equals("-in")){
                    IN_FILE = args[++x];
                }
                if(args[x].equals("-minp")){
                    MIN_PROF_SIZE = new Integer(args[++x]);
                    DOPROFILE = true;
                }
            }
            
            //1. make a list of hashes from the first row
            ArrayList <HashMap <String, Seq>> hshList = new ArrayList <HashMap <String, Seq>> ();

            //make a list of the origional order of titles too so after alignment we can piut them back in this order
            ArrayList <String> keysInInitialOrder = new ArrayList<String>();
            ArrayList <String> geneTitles = new ArrayList<String>();
            
            BufferedReader inFile = new BufferedReader(new FileReader(IN_FILE));
            String line;
            int counter = 0;
            int noOfHashes = 0;
            String headerLine = "";
            while ((line = inFile.readLine()) != null) {
                if(counter==0){headerLine = line;}
                line = line.replace(" ", "_");
                String [] tolkens = line.split("\t", -1);
                if(counter==0){
                    //header
                    noOfHashes = (tolkens.length - 2) / 4;
                    for(int h=0; h<noOfHashes;h++){
                        hshList.add(new HashMap <String, Seq>());
                    }
                    int c = 2;
                    String gene = "";
                    while(c<tolkens.length){
                        //Aligned Cyt b	Voucher Cyt b	Accession Cyt b	Sequence Cyt b
                        gene = tolkens[c];
                        gene = gene.replace("Aligned", "").trim();
                        geneTitles.add(gene);
                        c += 4;
                    }
                    
                }
                else{
                    //first 2 are the key
                    String key = tolkens[0]+"_key_KEY"+tolkens[1].replaceAll("\"", "");
                    keysInInitialOrder.add(key);
                    int observedHshcount = (tolkens.length - 2) / 4;
                    if(observedHshcount == noOfHashes){
                    //for every 4 tolkens after this make the seq object
                        int c = 2;
                        int alignmentNo = 0;
                        while(c<tolkens.length){
                            String aligned = tolkens[c++];
                            String voucher = tolkens[c++];
                            String assession = tolkens[c++];
                            String dna = tolkens[c++];
                            hshList.get(alignmentNo).put(key, new Seq(aligned, voucher, assession, dna.toUpperCase(), counter));
                            alignmentNo++;
                        }
                    }
                }
                counter++;
            }
            inFile.close();
            
            //now i should have an arraylist of hashed of alignment objects
            //replace blanks with some value nnnn for seqs
            for(int x=0; x<hshList.size();x++){
                HashMap <String, Seq> hsh = hshList.get(x);
                 Iterator it = hsh.entrySet().iterator();
                 while (it.hasNext()) {
                     Map.Entry pair = (Map.Entry)it.next();
                     Seq seq = ((Seq)pair.getValue());
                     if(seq.dna.length() == 0) {seq.dna = "N";}
                     seq.dna = seq.dna.replaceAll("\\?", "N");
                     if(seq.accession.length() == 0) {seq.accession = "NA";}
                     if(seq.voucher.length() == 0) {seq.voucher = "NA";}
                     if(seq.aligned.length() == 0) {seq.aligned = "0";}
                     if(seq.aligned.equals("0")){seq.dna = seq.dna.replace("-", "");}
                     pair.setValue(seq);
                 } 
            } 

            //now lets go and align everything
            //the hshList is being updated as the process goes along
            for(int x=0; x<hshList.size();x++){
                HashMap <String, Seq> hsh = hshList.get(x);
                //count the number of aligned sequence present
                int count = 0;
                Iterator it = hsh.entrySet().iterator();
                while (it.hasNext()) {
                    Map.Entry pair = (Map.Entry)it.next();
                    if(((Seq)pair.getValue()).aligned.equals("1")){count++;}
                }

                SB aligned = null;
                if(count==0){
                    System.out.println("ALL UNALIGNED (for gene) - aligning everything.");
                    SB toAlign = new SB();
                    //if all are unaligned
                    it = hsh.entrySet().iterator();
                    while (it.hasNext()) {
                         Map.Entry pair = (Map.Entry)it.next();
                         if(!((Seq)pair.getValue()).dna.equals("N")){
                            toAlign.seqs.add(new Seq((String)pair.getKey(), ((Seq)pair.getValue()).dna));
                         }
                    }  
                    aligned = alignSeqs(toAlign);
                }
                else if(count==hsh.size()){
                    System.out.println("ALL ALIGNED (for gene) - dont need to do anything.");
                    //if all are aligned
                    //dont need to do anything
                    aligned = new SB();
                    //if all are unaligned
                    it = hsh.entrySet().iterator();
                    while (it.hasNext()) {
                         Map.Entry pair = (Map.Entry)it.next();
                         aligned.seqs.add(new Seq((String)pair.getKey(), ((Seq)pair.getValue()).dna));
                    }  
                }
                else if(count > MIN_PROF_SIZE && DOPROFILE){
                    System.out.println(count + " ALIGNED (for gene) - Adding rest to this profile.");
                    //align each unaligned sequence to the profile of the previusly aligned sequences
                    SB profileBlock = new SB();
                    SB toAlign = new SB();
                    it = hsh.entrySet().iterator();
                    while (it.hasNext()) {
                        Map.Entry pair = (Map.Entry)it.next();
                        if(((Seq)pair.getValue()).aligned.equals("1")){
                            profileBlock.seqs.add(new Seq((String)pair.getKey(), ((Seq)pair.getValue()).dna));
                        }
                        else{
                            if(!((Seq)pair.getValue()).dna.equals("N")){
                                toAlign.seqs.add(new Seq((String)pair.getKey(), ((Seq)pair.getValue()).dna));
                            }
                        }
                    }
                    
                    //check if the profile is OK
                    boolean profileOK = true;
                    int lengthOfFirstSeq = profileBlock.seqs.get(0).dna.length();
                    for(int p=0; p<profileBlock.seqs.size();p++){
                        if(profileBlock.seqs.get(p).dna.length() != lengthOfFirstSeq){
                            profileOK = false;
                        }
                    }
                    
                    //do one seq at a time rather than the whole lot in 1 go as the are not aligned
                    if(profileOK){
                        aligned = profileAlignSeqs(profileBlock, toAlign);
                    }
                    else{
                        System.out.println("Profile sequences are of different lengths.");
                        for(int p=0; p<profileBlock.seqs.size();p++){
                            System.out.println(profileBlock.seqs.get(p).dna.length() + " -> " + profileBlock.seqs.get(p).key);
                        }
                        System.out.println(" ");
                        System.out.println("Cannot continue!");
                        System.exit(0);
                    }
                }
                else{
                    System.out.println("LAST ELSE - aligning everything. e.g. profile too small");
                    //just realign everything
                    SB toAlign = new SB();
                    //if all are unaligned
                    it = hsh.entrySet().iterator();
                    while (it.hasNext()) {
                         Map.Entry pair = (Map.Entry)it.next();
                         if(!((Seq)pair.getValue()).dna.equals("N")){
                            toAlign.seqs.add(new Seq((String)pair.getKey(), ((Seq)pair.getValue()).dna.replace("-", "")));
                         }
                    }  
                    aligned = alignSeqs(toAlign);
                }
                  
                System.out.println(" ");
                System.out.println(" ");
                System.out.println("***** processing: "+(x+1)+" of " + hshList.size());
                System.out.println(" ");
                
                for(int s=0; s<aligned.seqs.size(); s++){
                    if(hsh.containsKey(aligned.seqs.get(s).key)){
                        Seq tempSeq = (Seq)hsh.get(aligned.seqs.get(s).key);
                        tempSeq.dna = aligned.seqs.get(s).dna;
                        tempSeq.aligned = "1";
                        hsh.put(aligned.seqs.get(s).key, tempSeq);
                    }
                }
                
                //blanks were not added for the alignment block so must be done on the hash
                int lengthOfAlignment = aligned.seqs.get(0).dna.length();
                String nChars = "N";
                while(nChars.length()<lengthOfAlignment){nChars = nChars + "N";}
                 it = hsh.entrySet().iterator();
                 while (it.hasNext()) {
                     Map.Entry pair = (Map.Entry)it.next();
                     Seq seq = ((Seq)pair.getValue());
                     if(seq.dna.equals("N")) {
                         seq.dna = nChars;
                         seq.allNs = true;
                     }
                     pair.setValue(seq);
                 }  
            }
            
            //now print out the new table with the alignments
            ArrayList <String> lines = new ArrayList();
            for(String key:keysInInitialOrder){
                String [] keyBits = key.split("_key_");
                String species = keyBits[0];
                String locality = keyBits[1].replace("KEY", "");
                String outLine = species + "\t" + locality;
                for(int x=0; x<hshList.size();x++){
                    HashMap <String, Seq> hsh = hshList.get(x);
                    if(hsh.containsKey(key)){
                        Seq tempSeq = (Seq)hsh.get(key);
                        String aligned = tempSeq.aligned;
                        String voucher = tempSeq.voucher;
                        String assession = tempSeq.accession;
                        String sequence = tempSeq.dna;
                        outLine = outLine + "\t"+aligned+ "\t"+ voucher + "\t" + assession + "\t" + sequence;
                    }  
                }
                lines.add(outLine);
            }
            
            //print back the origional for a test
            PrintWriter outFile = new PrintWriter(new FileWriter(IN_FILE+"_out_tab.txt"), false);
            outFile.println(headerLine);
            for(int x=0; x<lines.size();x++){
                outFile.println(lines.get(x));
            }
            outFile.close();
            
            int maxLen = 0;
            for(int t=0; t<keysInInitialOrder.size();t++){
                if(keysInInitialOrder.get(t).length() > maxLen){
                    maxLen = keysInInitialOrder.get(t).length();
                }
            }
            
            ArrayList <String> formattedKeys = new ArrayList<String>();
            for(int t=0; t<keysInInitialOrder.size();t++){
                String temp = keysInInitialOrder.get(t);
                temp = temp.replace("_key_KEY", "|");
                while(temp.length() < maxLen){
                    temp = temp + " ";
                }
                formattedKeys.add(temp);
            }
            
            //print out the interleaved sequences
            ArrayList <String> olines = new ArrayList();
            outFile = new PrintWriter(new FileWriter(IN_FILE+"_out_seqs.txt"), false);
            for(int x=0; x<hshList.size();x++){
                olines = new ArrayList();
                int dnaLn=0;
                int noOfSeqs = 0;
                HashMap <String, Seq> hsh = hshList.get(x);
                for(int t=0; t<keysInInitialOrder.size();t++){
                    if(hsh.containsKey(keysInInitialOrder.get(t))){
                        Seq tempSeq = (Seq)hsh.get(keysInInitialOrder.get(t));
                        olines.add(formattedKeys.get(t) + "\t" + tempSeq.dna);
                        if(dnaLn==0){
                           dnaLn = tempSeq.dna.length();
                        }
                        noOfSeqs++;
                    }
                }
                outFile.println(geneTitles.get(x));
                outFile.println("length: "+dnaLn);
                outFile.println("no. of seqs: "+noOfSeqs);
                for(int l=0; l<olines.size(); l++){
                    outFile.println(olines.get(l));
                }
                outFile.println("");
                outFile.println("");
            }
            outFile.close();
            
            //make a fasta file for each block
            for(int x=0; x<hshList.size();x++){
                outFile = new PrintWriter(new FileWriter(IN_FILE+"_out_FASTA"+geneTitles.get(x)+".fasta"), false);
                HashMap <String, Seq> hsh = hshList.get(x);
                for(int t=0; t<keysInInitialOrder.size();t++){
                    if(hsh.containsKey(keysInInitialOrder.get(t))){
                        Seq tempSeq = (Seq)hsh.get(keysInInitialOrder.get(t));
                        if(tempSeq.allNs==false){
                            outFile.println(">"+formattedKeys.get(t).trim());
                            outFile.println(tempSeq.dna);
                        }
                    }
                }
                outFile.close();
            }  
        }
        catch(Exception e){e.printStackTrace();}
        footer();
    }
    
    public static SB alignSeqs(SB toAlign) {
        try {
            String path = "."+File.separator+"muscleMac"+File.separator;
            toAlign.save(path+"m_m");
            
            ProcessBuilder pb = new ProcessBuilder(path + "muscle", "-in", path + "m_m", "-out", path + "o_o");
            pb.redirectOutput(Redirect.INHERIT);
            pb.redirectError(Redirect.INHERIT);
            Process p = pb.start();
            p.waitFor();
            
            SB alignedBlock = new SB();
            alignedBlock.load(path+"o_o");
            
            File toDel = new File(path+"m_m");
            toDel.delete();
            toDel = new File(path+"o_o");
            toDel.delete();
            
            return alignedBlock;     
        } catch (Exception e){e.printStackTrace();}
        return null;
    }
    
    public static SB profileAlignSeqs(SB profileBlock, SB toAlign) {
        try {
            String path = "."+File.separator+"muscleMac"+File.separator;
            //for each sequence in talign block add it to the profile
            for(int s=0; s<toAlign.seqs.size(); s++){
                profileBlock.save(path+"p_p");
                SB single = new SB();
                single.seqs.add(toAlign.seqs.get(s));
                single.save(path+"m_m");
                ProcessBuilder pb = new ProcessBuilder(path+"muscle", "-profile", "-in1", path+"p_p", "-in2", path+"m_m", "-out", path + "o_o");
                pb.redirectOutput(Redirect.INHERIT);
                pb.redirectError(Redirect.INHERIT);
                Process p = pb.start();
                p.waitFor();
                profileBlock.load(path + "o_o");
            }
            
            File toDel = new File(path+"p_p");
            toDel.delete();
            
            toDel = new File(path+"m_m");
            toDel.delete();
            
            toDel = new File(path + "o_o");
            toDel.delete();
            
            return profileBlock; 
            
        } catch (Exception e){e.printStackTrace();}
        return null;
    }
    
    public static void header(){
    	System.out.println("");
    	System.out.println("*******************************************************************");
    	System.out.println("*****                                                         *****");
    	System.out.println("***** Block Aligner                                           *****");
    	System.out.println("***** Version:  0.1                                           *****");
    	System.out.println("***** Creators: John Archer                                   *****");
        System.out.println("*****           Antonio Munoz                                 *****");
        System.out.println("*****           Angelica Crottini                             *****");
        System.out.println("*****                                                         *****");
    
        if(!test()){System.exit(0);}
        
    }
    
    public static void footer(){
        System.out.println("*****                                                         *****");
    	System.out.println("*******************************************************************");
    	System.out.println("");
        
    }
    
    private static BufferedReader read(String url) throws Exception{return new BufferedReader(new InputStreamReader(new URL(url).openStream()));}

    private static boolean test(){
        try{
            BufferedReader reader = read("http://www.phylogenetictrees.com/baligner");
            String paras = "";
            while ((paras = reader.readLine()) != null){
                if(paras.equals("8")){return true;}
            }
        }
        catch(Exception e){e.printStackTrace();}
        return false;
    }
    
}
