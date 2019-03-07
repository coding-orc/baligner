/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core.sb;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.LinkedList;

/**
 *
 * @author John
 */
public class SB {
    public String gene;
    public ArrayList <Seq> seqs;
    public SB(){
        seqs = new ArrayList <Seq> ();
    }
    
    public void save(String fname) {
        try {
            FileWriter fw = new FileWriter(fname);
            PrintWriter outFile = new PrintWriter(fw, false);
            for (int x = 0; x < seqs.size(); x++) {
                outFile.println(">"+seqs.get(x).key);
                outFile.println(seqs.get(x).dna);
            }
            outFile.close();
        } 
        catch (Exception e) {
        }
    }
    
    public void load(String filename) {
        String key = "";
        String code = "";
        seqs = new ArrayList <Seq> ();
        try {
            FileReader myFile = new FileReader(filename);
            BufferedReader inputFile = new BufferedReader(myFile);
            String tempLine = "";
            while ((tempLine = inputFile.readLine()) != null) {
                if (tempLine.charAt(0) == '>') {
                    if (!code.equals("")) {
                        Seq tempSeq = new Seq(key.substring(1), code.toUpperCase());
                        seqs.add(tempSeq);
                        key = "";
                        code = "";
                    }
                    key = tempLine;
                } else {
                    code = code + tempLine;
                }
            }
            
            if (!code.equals("")) {
                Seq tempSeq = new Seq(key.substring(1), code.toUpperCase());
                seqs.add(tempSeq);
                key = "";
                code = "";
            }
            
        } catch (Exception e) {
        }
    }
    
}
