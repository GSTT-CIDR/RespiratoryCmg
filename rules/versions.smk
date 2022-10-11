rule getVersions:
    output:
        txt = "results/{sample}/{time}_hours/versions.txt"
    shell:
        """
        minimap2 --version | echo "minimap2 $(cat -)" >> {output.txt}
        centrifuge --version | head -n 1 | sed "s/.*bin\///g" - >> {output.txt}
        abricate --version >> {output.txt}
        scagaire --verison >> {output.txt}
        """

