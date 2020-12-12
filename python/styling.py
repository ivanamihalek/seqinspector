

""" Well ... styling

Part of seqinspector toy NGS QC pipeline
Ivana Mihalek,  2020

"""

class Styling:

	def __init__ (self, workbook=None):

		self.xlsx_format = {}
		if workbook:
			self.xlsx_format = {"header":workbook.add_format({'align': 'center',  'valign': 'vcenter',
															  'bold': True, 'text_wrap': True}),
								"red_border":workbook.add_format({'border': 5, 'border_color': 'red'}),
								"wordwrap":workbook.add_format({'align': 'left', 'text_wrap': True}),
								"hyperlink":workbook.add_format({'align': 'center', 'color': 'blue',
																 'underline': 1, 'valign': 'vcenter'})}
