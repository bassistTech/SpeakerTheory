'''
Crude markdown-to-HTML converter.

Francis Deck
MIT License
'''

import markdown
import sys

preamble = '''
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Max Width Example</title>
    <style>
        .container {
            max-width: 900px; /* Sets the maximum width of the container */
            width: 100%;      /* Ensures it is responsive and can be smaller than 900px */
            margin: 0 auto;   /* Centers the container on the page */
            background-color: #f0f0f0;
            padding: 20px;
        }
    </style>
</head>
<body>
  <div class="container">
'''

postamble = '''
  </div>
</body>
'''

def process_md_to_html(infile, outfile):
    md = open(infile).read()
    html_body = markdown.markdown(md)
    open(outfile, 'w').write(preamble + html_body + postamble)

if len(sys.argv) != 3:
    print('Need 2 arguments, input and output file')
else:
    infile = sys.argv[1]
    outfile = sys.argv[2]
    process_md_to_html(infile, outfile)