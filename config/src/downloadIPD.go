package main

import (
	"bufio"
	"bytes"
	"encoding/json"
	"fmt"
	"io"
	"net/http"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"sync"
	"time"

	"github.com/gofrs/flock"
)

func determineStartingPoint(dir string, gene string, date string) (int, []int) {

	var allValues []int

	// Define the path to the JSON file
	jsonFilePath := dir + "/" + gene + "_date_lookup.json"

	// Check if the file exists
	if _, err := os.Stat(jsonFilePath); os.IsNotExist(err) {
		return 1, allValues
	}

	// Read the JSON file
	jsonFile, err := os.Open(jsonFilePath)
	if err != nil {
		fmt.Println("Error opening JSON file:", err)
		return 1, allValues
	}
	defer jsonFile.Close()

	// Parse the JSON into a map
	byteValue, _ := io.ReadAll(jsonFile)
	var dateLookup map[string]int
	json.Unmarshal(byteValue, &dateLookup)

	// Collect all integer values in a slice
	for _, value := range dateLookup {
		allValues = append(allValues, value)
	}

	// Look for the date in the map
	if val, ok := dateLookup[date]; ok {
		return val, allValues
	}

	// If the date was not found
	return 1, allValues
}

func buildDateLookup(id string, dateStr string, lookupDir string) error {

	trimmedID := id[3:]
	intID, err := strconv.Atoi(trimmedID)
	if err != nil {
		return fmt.Errorf("failed to convert identifier to integer: %v", err)
	}

	jsonFilePath := filepath.Join(lookupDir, "date_lookup.json")

	// Initialize file lock
	fileLock := flock.New(jsonFilePath)

	// Try to lock with a retry mechanism
	retryCount := 0
	maxRetries := 5
	sleepDuration := 2 * time.Second

	var locked bool
	var lockErr error
	for retryCount < maxRetries {
		locked, lockErr = fileLock.TryLock()
		if locked {
			break
		}
		if lockErr != nil {
			return fmt.Errorf("failed to lock the file: %v", lockErr)
		}

		// Sleep for a while before retrying
		time.Sleep(sleepDuration)
		retryCount++
	}

	if !locked {
		return fmt.Errorf("could not obtain lock after %d retries", maxRetries)
	}

	defer fileLock.Unlock()

	dateMap := make(map[string]int)

	jsonFile, err := os.Open(jsonFilePath)
	if err != nil && !os.IsNotExist(err) {
		return fmt.Errorf("failed to open JSON file: %v", err)
	}
	defer jsonFile.Close()

	if jsonFile != nil {
		byteValue, _ := io.ReadAll(jsonFile)
		json.Unmarshal(byteValue, &dateMap)
	}

	dateMap[dateStr] = intID

	jsonData, err := json.MarshalIndent(dateMap, "", "  ")
	if err != nil {
		return fmt.Errorf("failed to marshal JSON: %v", err)
	}

	savePath := jsonFilePath
	if os.IsNotExist(err) {
		savePath = "date_lookup.json"
	}

	err = os.WriteFile(savePath, jsonData, 0644)
	if err != nil {
		return fmt.Errorf("failed to write JSON file: %v", err)
	}

	return nil
}

func defineUrls(startingPoint int, maxCount int, gene string) ([]string, []string) {

	sliceSize := maxCount - startingPoint + 1

	var urls = make([]string, sliceSize)
	var ids = make([]string, sliceSize)

	if gene == "MHC" || gene == "KIR" || gene == "HLA" || gene == "MHCPRO" || gene == "KIRPRO" {

		baseURL := ""
		prefix := "NHP"
		if gene == "MHC" {
			baseURL = "https://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=ipdmhc;id="
		} else if gene == "KIR" {
			baseURL = "https://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=ipdnhkir;id="
		} else if gene == "HLA" {
			baseURL = "https://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=imgthla;id="
			prefix = "HLA"
		} else if gene == "MHCPRO" {
			baseURL = "https://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=ipdmhcpro;id="
		} else if gene == "KIRPRO" {
			baseURL = "https://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=ipdnhkirpro;id="
		}

		for i := startingPoint; i <= maxCount; i++ {
			id := fmt.Sprintf("%s%05d", prefix, i)
			url := baseURL + id + ";style=raw"

			ids[i-startingPoint] = id
			urls[i-startingPoint] = url
		}
	}

	return urls, ids

}

// checkAndSaveFile reads the date 'DT' line of the given file, parses the date, and saves the file if the date is after dateStr.
func checkAndSaveFile(id string, lookupIDs []string, fileContent io.Reader, dateStr string, lookup_dir string) error {
	// Parse the given date string
	givenDate, err := time.Parse("2006-01-02", dateStr)
	if err != nil {
		return err
	}

	scanner := bufio.NewScanner(fileContent)

	// collect the idline and dateline
	var dateLine string
	var idLine string
	var fileContentBuffer bytes.Buffer
	for scanner.Scan() {
		line := scanner.Text()
		fileContentBuffer.WriteString(line + "\n") // Save the content to a buffer
		if strings.HasPrefix(line, "ERROR") {
			idLine = line
			dateLine = line
			break
		}
		if strings.HasPrefix(line, "ID") {
			// fmt.Println("Found line starting with 'ID':", line)
			idLine = line
		} else if strings.HasPrefix(line, ">") {
			idLine = line
		}
		if strings.HasPrefix(line, "DT") {
			// fmt.Println("Found line starting with 'DT':", line)
			dateLine = line
		}
	}

	if idLine == "ERROR 12 No entries found." {
		return nil
	}

	if strings.HasPrefix(idLine, ">") {

		tmp_str := strings.Split(idLine, " ")[0]
		idInFileStr := strings.Split(tmp_str, ":")[1]
		// define the output file name for a protein FASTA
		outputFilename := fmt.Sprintf("%s.fasta", idInFileStr)
		outputFile, err := os.Create(outputFilename)
		if err != nil {
			return err
		}
		defer outputFile.Close()
		_, err = io.Copy(outputFile, &fileContentBuffer)
		return err
	}

	// parse the date line
	dateLine = strings.TrimPrefix(dateLine, "DT")
	dateLine = strings.TrimSpace(dateLine)       // Remove any leading or trailing spaces
	dateInFileStr := strings.Fields(dateLine)[0] // Get the first field, which should be the date
	dateInFile, err := time.Parse("02/01/2006", dateInFileStr)
	if err != nil {
		return err
	}

	// parse the ID line
	idLine = strings.TrimPrefix(idLine, "ID")
	idLine = strings.TrimSpace(idLine)
	idInFileStr := strings.Split(idLine, ";")[0]
	idInFileStr = strings.TrimSpace(idInFileStr)

	// define the output file name
	outputFilename := fmt.Sprintf("%s.embl", idInFileStr)

	// if an id isn't in the date lookup, add to it
	found := false
	for _, str := range lookupIDs {
		if str == idInFileStr {
			found = true
			break
		}
	}
	if !found {
		err := buildDateLookup(idInFileStr, dateInFileStr, lookup_dir)
		if err != nil {
			return err
		}
	}

	// Compare the dates and save the file if the date in the file is after the given date
	if dateInFile.After(givenDate) {
		outputFile, err := os.Create(outputFilename)
		if err != nil {
			return err
		}
		defer outputFile.Close()
		_, err = io.Copy(outputFile, &fileContentBuffer)
		return err
	}

	return nil
}

func downloadFile(url string, id string, lookup_records []int, wg *sync.WaitGroup, dateStr string, lookup_dir string) error {
	// defer wg.Done() // Decrease counter when the goroutine completes

	var resp *http.Response
	var err error

	// Retry up to 3 times
	for retries := 0; retries < 3; retries++ {
		resp, err = http.Get(url)
		if err == nil && resp.StatusCode == http.StatusOK {
			break // Success, exit the retry loop
		}
		if err != nil {
			fmt.Printf("Error fetching URL: %s. Retrying...\n", err)
		} else {
			fmt.Printf("Received unexpected status code: %d. Retrying...\n", resp.StatusCode)
		}
		time.Sleep(2 * time.Second) // Wait before retrying
	}

	if err != nil {
		return fmt.Errorf("FAILED TO FETCH URL AFTER RETRIES: %w", err)
	}
	defer resp.Body.Close()

	if resp.StatusCode != http.StatusOK {
		return fmt.Errorf("RECEIVED UNEXPECTED STATUS CODE: %d", resp.StatusCode)
	}

	// Create a new slice to store the string IDs
	var allIDs []string
	if len(lookup_records) > 0 {

		// define ID prefix
		prefix := "NHP"
		if strings.Contains(url, "imgthla") {
			prefix = "HLA"
		}

		// Loop through lookup_records and convert each integer to a string ID
		for _, i := range lookup_records {
			id := fmt.Sprintf("%s%05d", prefix, i)
			allIDs = append(allIDs, id)
		}
	}

	// Check and save the file if conditions are met
	return checkAndSaveFile(id, allIDs, resp.Body, dateStr, lookup_dir)
}

func main() {

	// Check if enough arguments are provided
	if len(os.Args) < 4 {
		fmt.Println("Usage: download-ipd-alleles.go <gene> <number to download> <last IPD release date>")
		return
	}

	// Assign command-line arguments to variables
	gene := os.Args[1]
	alleleCount, err := strconv.Atoi(os.Args[2])
	if err != nil {
		fmt.Printf("Error converting '%s' to an integer: %v\n", os.Args[2], err)
		return
	}
	lastReleaseDate := os.Args[3]
	lookup_dir := os.Args[4]

	// Define a struct to hold URL and ID
	type UrlWithID struct {
		Url string
		ID  string
	}

	startingPoint, values := determineStartingPoint(lookup_dir, gene, lastReleaseDate)

	// retrieve IPD urls and ids
	urls, ids := defineUrls(startingPoint, alleleCount, gene)

	// Define the number of concurrent workers
	const numWorkers = 100

	// Create a channel to send UrlWithID structs to workers
	urlsChannel := make(chan UrlWithID, numWorkers)

	// Create a wait group to wait for all workers to complete
	var wg sync.WaitGroup

	// Launch workers
	for i := 0; i < numWorkers; i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for urlWithID := range urlsChannel {
				url, id := urlWithID.Url, urlWithID.ID
				err := downloadFile(url, id, values, &wg, lastReleaseDate, lookup_dir)
				if err != nil {
					fmt.Printf("Error downloading URL %s: %v\n", url, err)
				}
			}
		}()
	}

	// Send UrlWithID structs to the channel
	for i, url := range urls {
		urlWithID := UrlWithID{Url: url, ID: ids[i]}
		urlsChannel <- urlWithID
	}
	close(urlsChannel)

	// Wait for all workers to complete
	wg.Wait()

	fmt.Println("Download completed.")
}
