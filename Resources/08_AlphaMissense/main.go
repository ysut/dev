package main

import (
	"bufio"
	"compress/gzip"
	"database/sql"
	"fmt"
	"io"
	"log"
	"net/http"
	"os"
	"strings"

	_ "github.com/mattn/go-sqlite3"
)

// Mapping of table names to column names
// The key is the file name and the value is a list of table name and column names
var tableColumns = map[string][]string{
	"AlphaMissense_hg19.tsv.gz":                      {"hg19", "CHROM, POS, REF, ALT, genome, uniprot_id, transcript_id, protein_variant, am_pathogenicity, am_class"},
	"AlphaMissense_hg38.tsv.gz":                      {"hg38", "CHROM, POS, REF, ALT, genome, uniprot_id, transcript_id, protein_variant, am_pathogenicity, am_class"},
	"AlphaMissense_gene_hg19.tsv.gz":                 {"gene_hg19", "transcript_id, mean_am_pathogenicity"},
	"AlphaMissense_gene_hg38.tsv.gz":                 {"gene_hg38", "transcript_id, mean_am_pathogenicity"},
	"AlphaMissense_isoforms_hg38.tsv.gz":             {"isf_hg38", "CHROM, POS, REF, ALT, genome, transcript_id, protein_variant, am_pathogenicity, am_class"},
	"AlphaMissense_isoforms_aa_substitutions.tsv.gz": {"isf_aasub", "transcript_id, protein_variant, am_pathogenicity, am_class"},
	"AlphaMissense_aa_substitutions.tsv.gz":          {"aasub", "uniprot_id, protein_variant, am_pathogenicity, am_class"},
}

// Mapping of file names to URLs
var urls = map[string]string{
	"AlphaMissense_hg19.tsv.gz":                      "https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg19.tsv.gz",
	"AlphaMissense_hg38.tsv.gz":                      "https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz",
	"AlphaMissense_gene_hg19.tsv.gz":                 "https://storage.googleapis.com/dm_alphamissense/AlphaMissense_gene_hg19.tsv.gz",
	"AlphaMissense_gene_hg38.tsv.gz":                 "https://storage.googleapis.com/dm_alphamissense/AlphaMissense_gene_hg38.tsv.gz",
	"AlphaMissense_isoforms_hg38.tsv.gz":             "https://storage.googleapis.com/dm_alphamissense/AlphaMissense_isoforms_hg38.tsv.gz",
	"AlphaMissense_isoforms_aa_substitutions.tsv.gz": "https://storage.googleapis.com/dm_alphamissense/AlphaMissense_isoforms_aa_substitutions.tsv.gz",
	"AlphaMissense_aa_substitutions.tsv.gz":          "https://storage.googleapis.com/dm_alphamissense/AlphaMissense_aa_substitutions.tsv.gz",
}

// Check existence of a file
func Exists(filename string) bool {
	_, err := os.Stat(filename)
	return err == nil
}

// Download AlphaMissense bulk data
func DownloadFile(filepath string, url string) error {
	if Exists(filepath) {
		log.Println("File already exists: ", filepath)
		return nil
	} else {
		log.Println("Downloading file: ", filepath)
		// Get the data
		resp, err := http.Get(url)
		if err != nil {
			return err
		}
		defer resp.Body.Close()

		// Create the file
		out, err := os.Create(filepath)
		if err != nil {
			return err
		}
		defer out.Close()

		// Write the body to file
		_, err = io.Copy(out, resp.Body)
		return err
	}
}

func main() {
	// #0. Download the files
	for filePath, url := range urls {
		err := DownloadFile(filePath, url)
		if err != nil {
			log.Fatalf("Error downloading file: ", err)
		}
	}

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

	// #3. Start a transaction
	tx, err := db.Begin()
	if err != nil {
		log.Fatalf("Error starting transaction: %v", err)
	}
	defer tx.Rollback()

	// #4. Import data into the database defined in the dbPath
	for filePath, tableInfo := range tableColumns {
		importData(db, filePath, tableInfo[0], tableInfo[1])
	}

	// #5. Commit the transaction
	err = tx.Commit()
	if err != nil {
		log.Fatalf("Error committing transaction: %v", err)
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

	// Prepare the SQL statement for inserting data
	placeholders := strings.Repeat("?, ", len(strings.Split(columns, ", "))-1) + "?"
	insertSQL := fmt.Sprintf("INSERT INTO %s (%s) VALUES (%s)", tableName, columns, placeholders)
	stmt, err := db.Prepare(insertSQL)
	if err != nil {
		log.Fatalf("Error preparing statement: %v", err)
	}
	defer stmt.Close()

	// #3. Read the file line by line and insert into the database using prepared statement
	scanner := bufio.NewScanner(gz)

	// Skip first 4 lines (headers)
	for i := 0; i < 4; i++ {
		if !scanner.Scan() {
			log.Fatal("Failed to skip header lines")
		}
	}

	log.Println("Inserting data into", tableName)
	for scanner.Scan() {
		line := scanner.Text()
		values := strings.Split(line, "\t") // Split the line by tab

		// Convert slice of string to slice of interface{} for Exec
		interfaces := convertToInterfaceSlice(values)

		// Execute the prepared SQL statement
		_, err := stmt.Exec(interfaces...)
		if err != nil {
			log.Fatalf("Error inserting data into %s: %v", tableName, err)
		}
	}

	// Error checking for scanner
	if err := scanner.Err(); err != nil {
		log.Fatalf("Error reading from file %s: %v", filePath, err)
	}

	log.Printf("Data from %s imported into %s successfully.\n", filePath, tableName)
}

func convertToInterfaceSlice(slice []string) []interface{} {
	interfaces := make([]interface{}, len(slice))
	for i, v := range slice {
		interfaces[i] = v
	}
	return interfaces
}
