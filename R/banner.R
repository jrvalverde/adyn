#

charScanLine <- function(ch, l) {
    # This is the representation of each char as a series of 8 8-bit vectors
    #
    # To print a given character we must , for each of the eight lines
    # take the corresponding part and output the resulting horizontal scan 
    # line.
    #
    # This function will return, for a given character, the corresponding
    # horizontal scan line represented as a string of bits, where each
    # hash represents the corresponding bit
    #
    if      (ch == ' ')
        return(c(
        "        ", 
        "        ", 
        "        ", 
        "        ", 
        "        ", 
        "        ", 
        "        ", 
        "        "
	)[l])
    else if (ch == '0')
        return(c(
        " #####  ",
        "#    ## ",
        "#   # # ",
        "#  #  # ",
        "# #   # ",
        "##    # ",
        " #####  ",
        "        "
	)[l])
    else if (ch == "1")
        return(c(
        "     #  ",
        "    ##  ",
        "     #  ",
        "     #  ",
        "     #  ",
        "     #  ",
        "     #  ",
        "        "
        )[l])
    else if (ch == "2")
        return(c(
        " #####  ",
        "#     # ",
        "      # ",
        "  ####  ",
        " #      ",
        "#       ",
        "####### ",
        "        "
        )[l])
    else if (ch == "3")
        return(c(
        " #####  ",
        "#     # ",
        "      # ",
        "  ####  ",
        "      # ",
        "#     # ",
        " #####  ",
        "        "
        )[l])
    else if (ch == "4")
        return(c(
        "#    #  ",
        "#    #  ",
        "#    #  ",
        " ###### ",
        "     #  ",
        "     #  ",
        "     #  ",
        "        "
        )[l])
    else if (ch == "5")
        return(c(
        "####### ",
        "#       ",
        "#       ",
        "######  ",
        "      # ",
        "#     # ",
        " #####  ",
        "        "
        )[l])
    else if (ch == "6")
        return(c(
        " #####  ",
        "#       ",
        "#       ",
        " #####  ",
        "#     # ",
        "#     # ",
        " #####  ",
        "        "
        )[l])
    else if (ch == "7")
        return(c(
        " ###### ",
        "      # ",
        "     #  ",
        "    #   ",
        "   #    ",
        "  #     ",
        " #      ",
        "        "
        )[l])
    else if (ch == "8")
        return(c(
        " #####  ",
        "#     # ",
        "#     # ",
        " #####  ",
        "#     # ",
        "#     # ",
        " #####  ",
        "        "
        )[l])
    else if (ch == "9")
        return(c(
        " #####  ",
        "#     # ",
        "#     # ",
        " #####  ",
        "      # ",
        "      # ",
        " #####  ",
        "        "
        )[l])
    else if (ch == "a")
        return(c(
        "        ",
        "        ",
        " ####   ",
        "     #  ",
        " #####  ",
        "#    #  ",
        " #### # ",
        "        "      
	)[l])
    else if (ch == "b")
        return(c(
        " #      ",
        " #      ",
        " #####  ",
        " #    # ",
        " #    # ",
        " #    # ",
        "# ####  ",
        "        "
        )[l])
    else if (ch == "c")
        return(c(
        "        ",
        "        ",
        "  ###   ",
        " #   #  ",
        "#       ",
        " #   #  ",
        "  ###   ",
        "        "
        )[l])
    else if (ch == "d")
        return(c(
        "     #  ",
        "     #  ",
        " #####  ",
        "#    #  ",
        "#    #  ",
        "#    #  ",
        " #### # ",
        "        "
        )[l])
    else if (ch == "e")
        return(c(
        "        ",
        "        ",
        " ####   ",
        "#    #  ",
        "######  ",
        "#       ",
        " ####   ",
        "        "
        )[l])
    else if (ch == "f")
        return(c(
        "  ###   ",
        " #      ",
        " #      ",
        "####    ",
        " #      ",
        " #      ",
        " #      ",
        "        "        
        )[l])
    else if (ch == "g")
        return(c(
        "        ",
        "        ",
        " #### # ",
        "#    #  ",
        "#    #  ",
        " ####   ",
        "     #  ",
        " ####   "        
        )[l])
    else if (ch == "h")
        return(c(
        " #      ",
        " #      ",
        " ####   ",
        " #   #  ",
        " #   #  ",
        " #   #  ",
        "##   #  ",
        "        "
        )[l])
    else if (ch == "i")
        return(c(
        "  #     ",
        "        ",
        " ##     ",
        "  #     ",
        "  #     ",
        "  #     ",
        "#####   ",
        "        "
        )[l])
    else if (ch == "j")
        return(c(
        "   #    ",
        "        ",
        "  ##    ",
        "   #    ",
        "   #    ",
        "   #    ",
        "#  #    ",
        " ##     "        
        )[l])
    else if (ch == "k")
        return(c(
        " #      ",
        " #  #   ",
        " # #    ",
        " ##     ",
        " # #    ",
        " #  #   ",
        "##  ##  ",
        "        "
        )[l])
    else if (ch == "l")
        return(c(
        "##      ",
        " #      ",
        " #      ",
        " #      ",
        " #      ",
        " #      ",
        "###     ",
        "        "
        )[l])
    else if (ch == "m")
        return(c(
        "        ",
        "        ",
        "# # ##  ",
        "## #  # ",
        "#  #  # ",
        "#  #  # ",
        "#  #  # ",
        "        "
        )[l])
    else if (ch == "n")
        return(c(
        "        ",
        "        ",
        " # ##   ",
        " ##  #  ",
        " #   #  ",
        " #   #  ",
        " #   #  ",
        "        "
        )[l])
    else if (ch == "o")
        return(c(
        "        ",
        "        ",
        " ####   ",
        "#    #  ",
        "#    #  ",
        "#    #  ",
        " ####   ",
        "        "
        )[l])
    else if (ch == "p")
        return(c(
        "        ",
        "        ",
        "#####   ",
        " #   #  ",
        " #   #  ",
        " ####   ",
        " #      ",
        "##      "
        )[l])
    else if (ch == "q")
        return(c(
        "       ",
        "       ",
        " #### #",
        "#    # ",
        "#    # ",
        " ##### ",
        "     # ",
        "     ##" 
        )[l])
    else if (ch == "r")
        return(c(
        "        ",
        "        ",
        "## ##   ",
        " ##  #  ",
        " #      ",
        " #      ",
        "##      ",
        "        "
        )[l])
    else if (ch == "s")
        return(c(
        "        ",
        "        ",
        " ####   ",
        "#       ",
        " ####   ",
        "     #  ",
        "#####   ",
        "        "
        )[l])
    else if (ch == "t")
        return(c(
        "  #     ",
        "  #     ",
        "#####   ",
        "  #     ",
        "  #     ",
        "  #     ",
        "   ##   ",
        "        "
        )[l])
    else if (ch == "u")
        return(c(
        "        ",
        "        ",
        "#    #  ",
        "#    #  ",
        "#    #  ",
        "#    #  ",
        " #### # ",
        "        "        
        )[l])
    else if (ch == "v")
        return(c(
        "        ",
        "        ",
        "#     # ",
        " #   #  ",
        " #   #  ",
        "  # #   ",
        "   #    ",
        "        "        
        )[l])
    else if (ch == "w")
        return(c(
        "        ",
        "        ",
        "#     # ",
        "#     # ",
        "#  #  # ",
        "#  #  # ",
        " ## ##  ",
        "        "               
        )[l])
    else if (ch == "x")
        return(c(
        "        ",
        "        ",
        "#    #  ",
        " #  #   ",
        "  ##    ",
        " #  #   ",
        "#    #  ",
        "        "        
        )[l])
    else if (ch == "y")
        return(c(
        "        ",
        "        ",
        " #   #  ",
        " #   #  ",
        " #   #  ",
        "  ###   ",
        "    #   ",
        " ###    "        
        )[l])
    else if (ch == "z")
        return(c(
        "        ",
        "        ",
        "######  ",
        "    #   ",
        "   #    ",
        "  #     ",
        " ###### ",
        "        "        
        )[l])
    else if (ch == "A")
        return(c(
        "   #    ",
        "  # #   ",
        " #   #  ",
        "#     # ",
        "####### ",
        "#     # ",
        "#     # ",
        "        "
        )[l])
    else if (ch == "B")
        return(c(
        "#####   ",
        "#    #  ",
        "#    #  ",
        "#####   ",
        "#    #  ",
        "#    #  ",
        "#####   ",
        "        "
        )[l])
    else if (ch == "C")
        return(c(
        "  ###   ",
        " #   #  ",
        "#       ",
        "#       ",
        "#       ",
        " #   #  ",
        "  ###   ",
        "        "
        )[l])
    else if (ch == "D")
        return(c(
        "####    ",
        "#   #   ",
        "#    #  ",
        "#    #  ",
        "#    #  ",
        "#   #   ",
        "####    ",
        "        "
        )[l])
    else if (ch == "E")
        return(c(
        "######  ",
        "#       ",
        "#       ",
        "#####   ",
        "#       ",
        "#       ",
        "######  ",
        "        "
        )[l])
    else if (ch == "F")
        return(c(
        "######  ",
        "#       ",
        "#       ",
        "#####   ",
        "#       ",
        "#       ",
        "#       ",
        "        "
        )[l])
    else if (ch == "G")
        return(c(
        "  ###   ",
        " #   #  ",
        "#       ",
        "#       ",
        "#   ### ",
        " #   ## ",
        "  ### # ",
        "        "
        )[l])
    else if (ch == "H")
        return(c(
        "#     # ",
        "#     # ",
        "#     # ",
        "####### ",
        "#     # ",
        "#     # ",
        "#     # ",
        "        "
        )[l])
    else if (ch == "I")
        return(c(
        "#####   ",
        "  #     ",
        "  #     ",
        "  #     ",
        "  #     ",
        "  #     ",
        "#####   ",
        "        "
        )[l])
    else if (ch == "J")
        return(c(
        " #####  ",
        "   #    ",
        "   #    ",
        "   #    ",
        "   #    ",
        " # #    ",
        "  #     ",
        "        "
        )[l])
    else if (ch == "K")
        return(c(
        "#    #  ",
        "#   #   ",
        "#  #    ",
        "###     ",
        "#  #    ",
        "#   #   ",
        "#    #  ",
        "        "        
        )[l])
    else if (ch == "L")
        return(c(
        "#       ",
        "#       ",
        "#       ",
        "#       ",
        "#       ",
        "#       ",
        "######  ",
        "        "        
        )[l])
    else if (ch == "M")
        return(c(
        "##   ## ",
        "# # # # ",
        "#  #  # ",
        "#     # ",
        "#     # ",
        "#     # ",
        "#     # ",
        "        "        
        )[l])
    else if (ch == "N")
        return(c(
        "#     # ",
        "##    # ",
        "# #   # ",
        "#  #  # ",
        "#   # # ",
        "#    ## ",
        "#     # ",
        "        "        
        )[l])
    else if (ch == "O")
        return(c(
        "  ###   ",
        " #   #  ",
        "#     # ",
        "#     # ",
        "#     # ",
        " #   #  ",
        "  ###   ",
        "        "        
        )[l])
    else if (ch == "P")
        return(c(
        "#####   ",
        "#    #  ",
        "#    #  ",
        "#####   ",
        "#       ",
        "#       ",
        "#       ",
        "        "        
        )[l])
    else if (ch == "Q")
        return(c(
        "  ###   ",
        " #   #  ",
        "#     # ",
        "#     # ",
        "# #   # ",
        " # # #  ",
        "  ###   ",
        "     ## "        
        )[l])
    else if (ch == "R")
        return(c(
        "####    ",
        "#   #   ",
        "#   #   ",
        "####    ",
        "#   #   ",
        "#    #  ",
        "#    #  ",
        "        "        
        )[l])
    else if (ch == "S")
        return(c(
        "  ###   ",
        " #   #  ",
        "  #     ",
        "   #    ",
        "    #   ",
        " #   #  ",
        "  ###   ",
        "        "        
        )[l])
    else if (ch == "T")
        return(c(
        "####### ",
        "   #    ",
        "   #    ",
        "   #    ",
        "   #    ",
        "   #    ",
        "   #    ",
        "        "        
        )[l])
    else if (ch == "U")
        return(c(
        "#     # ",
        "#     # ",
        "#     # ",
        "#     # ",
        "#     # ",
        " #   #  ",
        "  ###   ",
        "        "
        )[l])
    else if (ch == "V")
        return(c(
        "#     # ",
        "#     # ",
        " #   #  ",
        " #   #  ",
        "  # #   ",
        "  # #   ",
        "   #    ",
        "        "
        )[l])
    else if (ch == "W")
        return(c(
        "#     # ",
        "#     # ",
        "#     # ",
        "#     # ",
        "#  #  # ",
        "# # # # ",
        " #   #  ",
        "        "
        )[l])
    else if (ch == "X")
        return(c(
        "#     # ",
        "#     # ",
        " #   #  ",
        "  ###   ",
        " #   #  ",
        "#     # ",
        "#     # ",
        "        "
        )[l])
    else if (ch == "Y")
        return(c(
        "#     # ",
        "#     # ",
        " #   #  ",
        "  # #   ",
        "   #    ",
        "   #    ",
        "   #    ",
        "        "
        )[l])
    else if (ch == "Z")
        return(c(
        "####### ",
        "     #  ",
        "    #   ",
        "   #    ",
        "  #     ",
        " #      ",
        "####### ",
        "        "
        )[l])
    else if (ch == ".")
        return(c(
        "        ",
        "        ",
        "        ",
        "        ",
        "        ",
        "        ",
        "   #    ",
        "        "        
        )[l])
    else if (ch == ",")
        return(c(
        "        ",
        "        ",
        "        ",
        "        ",
        "        ",
        "    #   ",
        "    #   ",
        "   #    "        
        )[l])
    else if (ch == ";")
        return(c(
        "        ",
        "        ",
        "        ",
        "    #   ",
        "        ",
        "    #   ",
        "    #   ",
        "   #    "        
        )[l])
    else if (ch == ":")
        return(c(
        "        ",
        "        ",
        "        ",
        "        ",
        "   #    ",
        "        ",
        "   #    ",
        "        "        
        )[l])
    else if (ch == "!")
        return(c(
        "   #    ",
        "   #    ",
        "   #    ",
        "   #    ",
        "   #    ",
        "        ",
        "   #    ",
        "        "        
        )[l])
    else if (ch == "\\")
        return(c(
        "#       ",
        " #      ",
        "  #     ",
        "   #    ",
        "    #   ",
        "     #  ",
        "      # ",
        "        "        
        )[l])
    else if (ch == "/")
        return(c(
        "      # ",
        "     #  ",
        "    #   ",
        "   #    ",
        "  #     ",
        " #      ",
        "#       ",
        "        "        
        )[l])
    else if (ch == "|")
        return(c(
        "   #    ",
        "   #    ",
        "   #    ",
        "   #    ",
        "   #    ",
        "   #    ",
        "   #    ",
        "   #    "        
        )[l])
    else if (ch == "+")
        return(c(
        "        ",
        "   #    ",
        "   #    ",
        " #####  ",
        "   #    ",
        "   #    ",
        "        ",
        "        "        
        )[l])
    else if (ch == "-")
        return(c(
        "        ",
        "        ",
        "        ",
        " #####  ",
        "        ",
        "        ",
        "        ",
        "        "        
        )[l])
    else if (ch == "*")
        return(c(
        "        ",
        " # # #  ",
        "  ###   ",
        "   #    ",
        "  ###   ",
        " # # #  ",
        "        ",
        "        "        
        )[l])
    else if (ch == "=")
        return(c(
        "        ",
        "        ",
        " #####  ",
        "        ",
        " #####  ",
        "        ",
        "        ",
        "        "        
        )[l])
    else if (ch == "_")
        return(c(
	"        ", 
        "        ", 
        "        ", 
        "        ", 
        "        ", 
        "        ", 
        "        ",
        "########"        
        )[l])
    else if (ch == "(")
        return(c(
	"     #  ", 
        "    #   ", 
        "   #    ", 
        "   #    ", 
        "   #    ", 
        "    #   ", 
        "     #  ",
        "        "       
        )[l])
    else if (ch == ")")
        return(c(
	"  #     ", 
        "   #    ", 
        "    #   ", 
        "    #   ", 
        "    #   ", 
        "   #    ", 
        "  #     ",
        "        "       
        )[l])
    else if (ch == "[")
        return(c(
	"   ###  ", 
        "   #    ", 
        "   #    ", 
        "   #    ", 
        "   #    ", 
        "   #    ", 
        "   ###  ",
        "        "       
        )[l])
    else if (ch == "]")
        return(c(
	"  ###   ", 
        "    #   ", 
        "    #   ", 
        "    #   ", 
        "    #   ", 
        "    #   ", 
        "  ###   ",
        "        "       
        )[l])
    else if (ch == "{")
        return(c(
	"    ##  ", 
        "   #    ", 
        "   #    ", 
        "  #     ", 
        "   #    ", 
        "   #    ", 
        "    ##  ",
        "        "       
        )[l])
    else if (ch == "}")
        return(c(
	"  ##    ", 
        "    #   ", 
        "    #   ", 
        "     #  ", 
        "    #   ", 
        "    #   ", 
        "  ##    ",
        "        "       
        )[l])
    else if (ch == "'")
        return(c(
        "    #   ",
        "    #   ",
        "   #    ",       
        "        ",
        "        ",
        "        ",
        "        ",
        "        "
        )[l])
    else if (ch == '"')
        return(c(
        "  # #   ",
        "  # #   ",
        "  # #   ",       
        "        ",
        "        ",
        "        ",
        "        ",
        "        "
        )[l])
    else               
        return(c(
        "########",
        "#      #",
        "#      #",
        "#      #",
        "#      #",
        "#      #",
        "#      #",
        "########"        
        )[l])

}

banner <- function(text=" ")
{
    #len <- nchar(text)
    chars <- strsplit(text, "")[[1]]
    
    for (l in 1:8) {
        for (ch in chars) {
            cat(charScanLine(ch, l))
        }
        cat('\n')
    }
}

#banner('_______')
#banner()
#banner(' Hiya! ')
#banner('_______')
