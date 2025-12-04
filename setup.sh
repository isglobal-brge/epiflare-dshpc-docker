#!/bin/bash

# Colors for output
GREEN="\033[1;32m"
BLUE="\033[1;34m"
YELLOW="\033[1;33m"
RED="\033[1;31m"
CYAN="\033[1;36m"
TURQUOISE="\033[0;36m"
ORANGE="\033[0;33m"
PURPLE="\033[0;35m"
BOLD="\033[1m"
NC="\033[0m" # No Color

# Configuration
SCRIPT_NAME="configure-environment.sh"
GITHUB_RAW_URL="https://raw.githubusercontent.com/isglobal-brge/dsHPC-core/main/${SCRIPT_NAME}"

# Print banner
echo -e "${BLUE}┌────────────────────────────────────────────────────┐${NC}"
echo -e "${BLUE}│${NC} ${BLUE}              __       __  __ ____   ______       ${NC} ${BLUE}│${NC}"
echo -e "${BLUE}│${NC} ${BLUE}         ____/ /_____ / / / // __ \ / ____/       ${NC} ${BLUE}│${NC}"
echo -e "${BLUE}│${NC} ${BLUE}        / __  // ___// /_/ // /_/ // /            ${NC} ${BLUE}│${NC}"
echo -e "${BLUE}│${NC} ${BLUE}       / /_/ /(__  )/ __  // ____// /____         ${NC} ${BLUE}│${NC}"
echo -e "${BLUE}│${NC} ${BLUE}       \__,_//____//_/ /_//_/     \_____/         ${NC} ${BLUE}│${NC}"
echo -e "${BLUE}│${NC} ${BLUE}                                                  ${NC} ${BLUE}│${NC}"
echo -e "${BLUE}└────────────────────────────────────────────────────┘${NC}"
echo -e ""
echo -e "${BOLD}Welcome to ${YELLOW}High-Performance Computing for DataSHIELD${NC}${BOLD}!${NC}"
echo -e "       Made with ❤️  by ${BOLD}\033]8;;https://davidsarratgonzalez.github.io\007David Sarrat González\033]8;;\007${NC}"
echo ""
echo -e "${TURQUOISE}${BOLD}\033]8;;https://brge.isglobal.org\007Bioinformatics Research Group in Epidemiology (BRGE)\033]8;;\007${NC}"
echo -e "  ${ORANGE}${BOLD}\033]8;;https://www.isglobal.org\007Barcelona Institute for Global Health (ISGlobal)\033]8;;\007${NC}"
echo ""
echo

echo -e "${CYAN}📥 Downloading latest configuration script...${NC}"
echo -e "   Source: ${GITHUB_RAW_URL}"
echo

# Try to download the configure script
download_status=1
if command -v curl &> /dev/null; then
    curl -fsSL "$GITHUB_RAW_URL" -o "$SCRIPT_NAME"
    download_status=$?
elif command -v wget &> /dev/null; then
    wget -q "$GITHUB_RAW_URL" -O "$SCRIPT_NAME"
    download_status=$?
else
    echo -e "${RED}❌ Error: Neither curl nor wget is available${NC}"
    echo -e "${YELLOW}Please install curl or wget to continue${NC}"
    exit 1
fi

# Check if download was successful
if [[ $download_status -ne 0 ]]; then
    # Download failed - check if local copy exists
    if [[ -f "$SCRIPT_NAME" ]]; then
        echo -e "${YELLOW}⚠️  Could not download from GitHub (not yet pushed or network issue)${NC}"
        echo -e "${GREEN}✓ Using existing local copy of ${SCRIPT_NAME}${NC}"
    else
        echo -e "${RED}❌ Failed to download ${SCRIPT_NAME} and no local copy found${NC}"
        echo -e "${YELLOW}Please check your internet connection or ensure the file is pushed to GitHub${NC}"
        exit 1
    fi
else
    echo -e "${GREEN}✓ Downloaded ${SCRIPT_NAME} from GitHub${NC}"
fi

echo

# Make the script executable
chmod +x "$SCRIPT_NAME"

# Execute the configuration script
echo -e "${CYAN}▶️  Executing configuration script...${NC}"
echo

./"$SCRIPT_NAME" "$@"
execution_status=$?

echo
if [[ $execution_status -eq 0 ]]; then
    echo -e "${GREEN}✅ Setup completed successfully!${NC}"
else
    echo -e "${RED}❌ Setup encountered errors${NC}"
    exit $execution_status
fi

echo
echo -e "${CYAN}💡 Note: ${SCRIPT_NAME} is automatically downloaded and can be safely deleted${NC}"
echo -e "${CYAN}   It will be re-downloaded on next setup run${NC}"
echo
