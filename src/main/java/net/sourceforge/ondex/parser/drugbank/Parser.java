package net.sourceforge.ondex.parser.drugbank;

import net.sourceforge.ondex.parser.ONDEXParser;
import net.sourceforge.ondex.args.ArgumentDefinition;
import net.sourceforge.ondex.args.FileArgumentDefinition;
import net.sourceforge.ondex.core.ONDEXConcept;
import net.sourceforge.ondex.core.ONDEXRelation;
import net.sourceforge.ondex.core.ConceptClass;
import net.sourceforge.ondex.core.RelationType;
import net.sourceforge.ondex.core.CV;
import net.sourceforge.ondex.core.EvidenceType;
import net.sourceforge.ondex.core.AttributeName;

import net.sourceforge.ondex.exception.type.MetaDataMissingException;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.Hashtable;
import java.util.Map;
import java.util.ArrayList;

public class Parser extends ONDEXParser {
    /* 
    Author - sjcockell
    Date - 20100310
    Version - 0.0.1
     */

    private ConceptClass drugCC;
    private RelationType proteinProteinRelation;
    private CV genericCV;
    private EvidenceType genericET;
    private AttributeName inchiGDS;
    private AttributeName formulaGDS;
    private AttributeName indicationGDS;
    private CV casCV;
    private CV smilesCV;
    private CV keggCV;
    private CV wikiCV;
    private CV inchiCV;
    private CV pubchemCV;
    private CV chebiCV;
    private ONDEXConcept MatchedConcept = null;
    private Map<String, ONDEXConcept> proteins = new Hashtable<String, ONDEXConcept>();
    private Map<String, String> uniprotInfo = new Hashtable<String, String>();
    private Map<String, String> hprdInfo = new Hashtable<String, String>();
    private Map<String, String> hugoInfo = new Hashtable<String, String>();
    private Map<String, String> gsInfo = new Hashtable<String, String>();
    private Map<String, String> gnInfo = new Hashtable<String, String>();
    private Map<String, String> affyMapping = new Hashtable<String, String>();
	private Map<String, Float> expressionData = new Hashtable<String, Float>();
    public boolean readsDirectory() {
        return false;
    }

    public boolean readsFile() {
        return true;
    }

    public ArgumentDefinition<?>[] getArgumentDefinitions() {
        /*
        new ArgumentDefinition<?>[] {
        new StringArgumentDefinition(params)
        }
         */
        return new ArgumentDefinition<?>[]{
            new FileArgumentDefinition(FileArgumentDefinition.INPUT_FILE, "File to import", true, true, true, false)
        };
    }

    public String getName() {
        return "DrugBank";
    }

    public String getId() {
        return "drugbank";
    }

    public String getVersion() {
        return "0.1.3";
    }

    public String[] requiresValidators() {
        return new String[0];
    }

    public boolean requiresIndexedGraph() {
        return false;
    }

    public void start() throws Exception {
        fetchMetadata();
        //get directory name
        String filename = (String)pa.getUniqueValue(FileArgumentDefinition.INPUT_FILE);
        //read primary file
        BufferedReader br = new BufferedReader(new FileReader(filename));
        //process the file, line at a time
        String line = null;
        int i = 0;
        String entry = null;
        while ((line = br.readLine()) != null) {
            //#BEGIN_DRUGCARD DB00001
            if (line.startsWith("#BEGIN_DRUGCARD")) {
                entry = line + "\n";
            }
            //#END_DRUGCARD DB00001
            else if (line.startsWith("#END_DRUGCARD")) {
                entry = entry + line + "\n";
                processDrugcard(entry, i);
                i += 1;
            }
            else {
                entry = entry + line + "\n";
            }
            /*    
            METHODS FOR ANNOTATING CONCEPTS
            createConceptAccession(String accession, CV, boolean ambiguous)
            createConceptName(String name, boolean preferred
            createConceptGDS(AttributeName attrname, Object value, boolean doIndex)
             */
            //c1 = graph.getFactory().createConcept(++i + "", genericCV, proteinCC, genericET);
            //c1.createConceptName(cID1, false);
            //c1.createConceptAccession(cID1, genericCV, false);
        }
    }

    private void processDrugcard(String e, int index) {
        ONDEXConcept c1 = null;
        String[] split = e.split("\n");
        ArrayList brands = new ArrayList();
        ArrayList interactions = new ArrayList();
        ArrayList acc_nos = new ArrayList();
        String generic_name = new String();
        String cas_number = new String();
        String chebi_id = new String();
        String formula = new String();
        String inchi_id = new String();
        String inchi_key = new String();
        String indication = new String();
        String kegg = new String();
        String pubchem = new String();
        String smiles = new String();
        String wikipedia = new String();
        for (int i = 0; i < split.length; i++) {
            String ln = split[i];
            if (ln.startsWith("#BEGIN_DRUGCARD")) {
                String db_id = ln.split(" ")[1].trim();
                c1 = graph.getFactory().createConcept(index+"", genericCV, drugCC, genericET);
                c1.createConceptAccession(db_id, genericCV, false);
            }
            else if (ln.startsWith("# ")) {
                String entry_name = ln.split(" ")[1].trim();
                if (entry_name.equals("Generic_Name:")) {
                    i += 1;
                    if (!split[i].trim().equals("Not Available")) {
                        generic_name = split[i].trim();
                        System.out.println(generic_name);
                        c1.createConceptName(generic_name, true);
                    }
                }
                if (entry_name.equals("Brand_Names:")) {
                    while (!ln.trim().equals("")) {
                        i += 1;
                        ln = split[i];
                        if (!ln.trim().equals("Not Available")) {
                            brands.add(ln.trim());
                        }
                    }
                    for (int brandint=0; brandint < brands.size(); brandint++) {
                        String brand_name = (String)brands.get(brandint);
                        if (brand_name.length() != 0) {
                            c1.createConceptName(brand_name, false);
                        }
                    }
                }
                else if (entry_name.equals("CAS_Registry_Number:")) {
                    i += 1;
                    if (!split[i].trim().equals("Not Available")) {
                        cas_number = split[i].trim();
                        c1.createConceptAccession(cas_number, casCV, false);
                    }
                }
                else if (entry_name.equals("ChEBI_ID:")) {
                    i += 1;
                    if (!split[i].trim().equals("Not Available")) {
                        chebi_id = split[i].trim();
                        c1.createConceptAccession(chebi_id, chebiCV, false);
                    }
                }
                else if (entry_name.equals("Chemical_Formula:")) {
                    i += 1;
                    if (!split[i].trim().equals("Not Available")) {
                        formula = split[i].trim();
                        c1.createGDS(formulaGDS, formula, false);
                        //System.out.println(formula);
                    }
                }
                else if (entry_name.equals("Drug_Interactions:")) {
                    while (!ln.trim().equals("")) {
                        i += 1;
                        ln = split[i];
                        if (!ln.trim().equals("Not Available")) {
                            interactions.add(ln.trim());
                        }
                    }
                }
                else if (entry_name.equals("InChI_Identifier:")) {
                    i += 1;
                    if (!split[i].trim().equals("Not Available")) {
                        inchi_id = split[i].trim();
                        c1.createGDS(inchiGDS, inchi_id, false);
                    }
                }
                else if (entry_name.equals("InChI_Key:")) {
                    i += 1;
                    if (!split[i].trim().equals("Not Available")) {
                        inchi_key = split[i].trim();
                        c1.createConceptAccession(inchi_key, inchiCV, false);
                    }
                }
                else if (entry_name.equals("Indication:")) {
                    i += 1;
                    if (!split[i].trim().equals("Not Available")) {
                        indication = split[i].trim();
                        c1.createGDS(indicationGDS, indication, false);
                        //System.out.println(indication);
                    }
                }
                else if (entry_name.equals("KEGG_Compound_ID:")) {
                    i += 1;
                    if (!split[i].trim().equals("Not Available")) {
                        kegg = split[i].trim();
                        c1.createConceptAccession(kegg, keggCV, false);
                    }
                }
                else if (entry_name.equals("PubChem_Compound_ID:")) {
                    i += 1;
                    if (!split[i].trim().equals("Not Available")) {
                        pubchem = split[i].trim();
                        c1.createConceptAccession(pubchem, pubchemCV, false);
                    }
                }
                else if (entry_name.equals("Secondary_Accession_No:")) {
                    while (!ln.trim().equals("")) {
                        i += 1;
                        ln = split[i];
                        if (!ln.trim().equals("Not Available")) {
                            acc_nos.add(ln.trim());
                        }
                    }
                }
                else if (entry_name.equals("Smiles_String_canonical:")) {
                    i += 1;
                    if (!split[i].trim().equals("Not Available")) {
                        smiles = split[i].trim();
                        c1.createConceptAccession(smiles, smilesCV, false);
                    }
                }
                else if (entry_name.equals("Wikipedia_Link:")) {
                    i += 1;
                    if (!split[i].trim().equals("Not Available")) {
                        wikipedia = split[i].trim();
                        c1.createConceptAccession(wikipedia, wikiCV, false);
                    }
                }

            }
            else {
            }
        }
    }

    private void fetchMetadata() throws MetaDataMissingException {
        drugCC = requireConceptClass("Comp");
        genericCV = requireCV("DRUGBANK");
        genericET = requireEvidenceType("IMPD");
        casCV = requireCV("CAS");
        smilesCV = requireCV("SMILES");
        inchiCV = requireCV("INCHI");
        keggCV = requireCV("KEGG");
        wikiCV = requireCV("WIKI");
        pubchemCV = requireCV("PUBCHEM");
        chebiCV = requireCV("CHEBI");
        inchiGDS = requireAttributeName("INCHI_STRING");
        formulaGDS = requireAttributeName("FORMULA");
        indicationGDS = requireAttributeName("INDICATION");
    }

}

