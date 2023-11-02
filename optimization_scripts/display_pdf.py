import os
from PyPDF2 import PdfMerger
import argparse

def merge_and_display_pdfs(pdf_directory, output_pdf_path):
    # Get a list of PDF files in the directory and sort them alphabetically
    pdf_files = sorted([file for file in os.listdir(pdf_directory) if file.endswith('.pdf')])

    # Merge PDF files into one
    merger = PdfMerger()
    for pdf_file in pdf_files:
        pdf_path = os.path.join(pdf_directory, pdf_file)
        merger.append(pdf_path)

    # Save the merged PDF to a file
    with open(output_pdf_path, 'wb') as output_pdf_file:
        merger.write(output_pdf_file)

    # Display the merged PDF file (requires a PDF viewer)
    os.system(f'xdg-open {output_pdf_path}')  # For Linux

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Merge and display all PDF files in a directory.')
    parser.add_argument('pdf_directory', type=str, help='Directory containing the PDF files')
    parser.add_argument('output_pdf_path', type=str, help='Output path for the merged PDF')
    args = parser.parse_args()

    merge_and_display_pdfs(args.pdf_directory, args.output_pdf_path)

