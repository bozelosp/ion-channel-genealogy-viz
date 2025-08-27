import { NextRequest, NextResponse } from 'next/server';
import * as Diff from 'diff';

function escapeHtml(text: string): string {
  const map: { [key: string]: string } = {
    '&': '&amp;',
    '<': '&lt;',
    '>': '&gt;',
    '"': '&quot;',
    "'": '&#039;'
  };
  return text.replace(/[&<>"']/g, (m) => map[m]);
}

function generateHtmlDiff(string1: string, string2: string): string {
  const lines1 = string1.split('\n');
  const lines2 = string2.split('\n');
  
  const diff = Diff.diffLines(string1, string2);
  
  let html = `
    <!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
    "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
    <html>
    <head>
      <meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
      <title></title>
      <style type="text/css">
        table.diff {font-family:Courier; border:medium;}
        .diff_header {background-color:#e0e0e0}
        td.diff_header {text-align:right}
        .diff_next {background-color:#c0c0c0}
        .diff_add {background-color:#aaffaa}
        .diff_chg {background-color:#ffff77}
        .diff_sub {background-color:#ffaaaa}
      </style>
    </head>
    <body>
    <table class="diff" id="difflib_chg_to0__top" cellspacing="0" cellpadding="0" rules="groups">
      <colgroup></colgroup> <colgroup></colgroup> <colgroup></colgroup>
      <colgroup></colgroup> <colgroup></colgroup> <colgroup></colgroup>
      <thead>
        <tr>
          <th class="diff_next"><br /></th>
          <th colspan="2" class="diff_header">Source Code 1</th>
          <th class="diff_next"><br /></th>
          <th colspan="2" class="diff_header">Source Code 2</th>
        </tr>
      </thead>
      <tbody>
  `;

  let line1 = 1;
  let line2 = 1;

  for (const part of diff) {
    const lines = part.value.split('\n');
    // Remove last empty line if it exists
    if (lines[lines.length - 1] === '') {
      lines.pop();
    }

    if (part.added) {
      // Added lines (green)
      for (const line of lines) {
        html += `
          <tr>
            <td class="diff_next" id="difflib_chg_to0__${line2}"></td>
            <td class="diff_header" id="from0_${line1}"></td>
            <td nowrap="nowrap"></td>
            <td class="diff_next"></td>
            <td class="diff_header" id="to0_${line2}">${line2}</td>
            <td nowrap="nowrap" class="diff_add">+&nbsp;${escapeHtml(line)}</td>
          </tr>
        `;
        line2++;
      }
    } else if (part.removed) {
      // Removed lines (red)
      for (const line of lines) {
        html += `
          <tr>
            <td class="diff_next" id="difflib_chg_to0__${line1}"></td>
            <td class="diff_header" id="from0_${line1}">${line1}</td>
            <td nowrap="nowrap" class="diff_sub">-&nbsp;${escapeHtml(line)}</td>
            <td class="diff_next"></td>
            <td class="diff_header" id="to0_${line2}"></td>
            <td nowrap="nowrap"></td>
          </tr>
        `;
        line1++;
      }
    } else {
      // Unchanged lines
      for (const line of lines) {
        html += `
          <tr>
            <td class="diff_next"></td>
            <td class="diff_header" id="from0_${line1}">${line1}</td>
            <td nowrap="nowrap">&nbsp;${escapeHtml(line)}</td>
            <td class="diff_next"></td>
            <td class="diff_header" id="to0_${line2}">${line2}</td>
            <td nowrap="nowrap">&nbsp;${escapeHtml(line)}</td>
          </tr>
        `;
        line1++;
        line2++;
      }
    }
  }

  html += `
      </tbody>
    </table>
    <table class="diff" summary="Legends">
      <tr> <th colspan="2"> Legends </th> </tr>
      <tr> <td> <table border="" summary="Colors">
                    <tr><th> Colors </th> </tr>
                    <tr><td class="diff_add">&nbsp;Added&nbsp;</td></tr>
                    <tr><td class="diff_chg">Changed</td> </tr>
                    <tr><td class="diff_sub">Deleted</td> </tr>
                  </table></td>
           <td> <table border="" summary="Links">
                    <tr><th colspan="2"> Links </th> </tr>
                    <tr><td>(f)irst change</td> </tr>
                    <tr><td>(n)ext change</td> </tr>
                    <tr><td>(t)op</td> </tr>
                  </table></td> </tr>
    </table>
    </body>
    </html>
  `;

  return html;
}

export async function GET(request: NextRequest) {
  try {
    const { searchParams } = new URL(request.url);
    const string1 = searchParams.get('string1');
    const string2 = searchParams.get('string2');

    if (!string1 || !string2) {
      return NextResponse.json(
        { error: 'Both string1 and string2 parameters are required' },
        { status: 400 }
      );
    }

    // Generate HTML diff
    const diffHtml = generateHtmlDiff(
      decodeURIComponent(string1),
      decodeURIComponent(string2)
    );

    return new NextResponse(diffHtml, {
      status: 200,
      headers: {
        'Content-Type': 'text/html',
      },
    });

  } catch (error) {
    console.error('Error generating diff:', error);
    return NextResponse.json(
      { error: `An error occurred: ${error instanceof Error ? error.message : 'Unknown error'}` },
      { status: 500 }
    );
  }
}