package main

import (
	"bufio"
	"compress/gzip"
	"database/sql"
	"fmt"
	"log"
	"os"
	"strings"

	_ "github.com/mattn/go-sqlite3"
)

// Mapping of table names to column names
// The key is the file name and the value is a list of table name and column names
var tableColumns = map[string][]string{
	"AlphaMissense_gene_hg19.tsv.gz":                 {"gene_hg19", "transcript_id, mean_am_pathogenicity"},
	"AlphaMissense_hg19.tsv.gz":                      {"hg19", "CHROM, POS, REF, ALT, genome, uniprot_id, transcript_id, protein_variant, am_pathogenicity, am_class"},
	"AlphaMissense_hg38.tsv.gz":                      {"hg38", "CHROM, POS, REF, ALT, genome, uniprot_id, transcript_id, protein_variant, am_pathogenicity, am_class"},
	"AlphaMissense_gene_hg38.tsv.gz":                 {"gene_hg38", "transcript_id, mean_am_pathogenicity"},
	"AlphaMissense_isoforms_hg38.tsv.gz":             {"isf_hg38", "CHROM, POS, REF, ALT, genome, transcript_id, protein_variant, am_pathogenicity, am_class"},
	"AlphaMissense_isoforms_aa_substitutions.tsv.gz": {"isf_aasub", "transcript_id, protein_variant, am_pathogenicity, am_class"},
	"AlphaMissense_aa_substitutions.tsv.gz":          {"aasub", "uniprot_id, protein_variant, am_pathogenicity, am_class"},
}

func main() {
	// #1. Open the database
	dbPath := "am_database.db" // ToDo: make it configurable using arguments
	db, err := sql.Open("sqlite3", dbPath)
	if err != nil {
		log.Fatalf("Error opening database: %v", err)
	}
	defer db.Close()

	// #2. Create tables in the database defined in the dbPath
	for _, info := range tableColumns {
		tableName := info[0]
		columns := info[1]
		createTableSQL := fmt.Sprintf("CREATE TABLE IF NOT EXISTS %s (%s)", tableName, columns)
		_, err := db.Exec(createTableSQL)
		if err != nil {
			log.Fatalf("Failed to create table %s: %v", tableName, err)
		}
	}

	// #3. Import data into the database defined in the dbPath
	for filePath, tableInfo := range tableColumns {
		importData(db, filePath, tableInfo[0], tableInfo[1])
	}
}

// Function to import data from a gzip file into a table
func importData(db *sql.DB, filePath, tableName, columns string) {
	// #1. Open the gzip file
	f, err := os.Open(filePath)
	if err != nil {
		log.Fatalf("Error opening file %s: %v", filePath, err)
	}
	defer f.Close()

	// #2. Create a gzip reader
	gz, err := gzip.NewReader(f)
	if err != nil {
		log.Fatalf("Error creating gzip reader for file %s: %v", filePath, err)
	}
	defer gz.Close()

	// #3. Read the file line by line and insert into the database
	scanner := bufio.NewScanner(gz)

	if scanner.Scan() {
		// pass the header
	}

	for scanner.Scan() {
		line := scanner.Text()
		values := strings.Split(line, "\t") // Split the line by tab

		// Generate the placeholders for the SQL query
		// If final "?" is not added,
		//     placeholder will be "?, ?, ... ?, " instead of "?, ?, ... ?"
		placeholders := strings.Repeat("?, ", len(values)-1) + "?"

		// Generate the SQL query
		insertSQL := fmt.Sprintf(
			"INSERT INTO %s (%s) VALUES (%s)", tableName, columns, placeholders)

		// Execute the SQL query
		// If there is an error,
		//     Fatalf will print the error and exit the program
		_, err := db.Exec(insertSQL, convertToInterfaceSlice(values)...)
		if err != nil {
			log.Fatalf("Error inserting data into %s: %v", tableName, err)
		}
	}

	// Error checking for scanner
	if err := scanner.Err(); err != nil {
		log.Fatalf("Error reading from file %s: %v", filePath, err)
	}

	fmt.Printf("Data from %s imported into %s successfully.\n", filePath, tableName)
}

func convertToInterfaceSlice(slice []string) []interface{} {
	interfaces := make([]interface{}, len(slice))
	for i, v := range slice {
		interfaces[i] = v
	}
	return interfaces
}
