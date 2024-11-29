from markdown_it import MarkdownIt

def convert_md_to_html(input_file, output_file):
    # Read the Markdown content
    with open(input_file, "r") as md_file:  # 'r' to read
        markdown_content = md_file.read()

    # Create a Markdown parser
    md = MarkdownIt()
    html_content = md.render(markdown_content)

    # Write the HTML content to an output file
    with open(output_file, "w") as html_file:   # 'w' to write
        html_file.write(html_content)

    print(f"HTML file created: {output_file}")  # Prints the successful conversion

# Converts the README from Markdown to HTML
convert_md_to_html("README.md", "README.html")
