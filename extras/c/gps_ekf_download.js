function offerFileAsDownload(filename, mime) {
	mime = mime || "application/octet-stream";

	let content = Module.FS.readFile(filename);
	console.log(`Offering download of "${filename}", with ${content.length} bytes...`);

	var a = document.createElement('a');
	a.download = filename;
	a.href = URL.createObjectURL(new Blob([content], {type: mime}));
	a.style.display = 'none';

	document.body.appendChild(a);
	a.click();
	setTimeout(() => {
		document.body.removeChild(a);
		URL.revokeObjectURL(a.href);
	}, 2000);
};

//roda e já entra baixando a saída.
if(Module){
	Module.onRuntimeInitialized = () => {
		Module.postRun.push(() => offerFileAsDownload("ekf.csv", "text/csv"));
	};
}
