# VoLTEVoWiFiHandoverPolicy
This is a simulation that discusses handover between VoLTE and VoWiFi

Simulation engine: Lab117
	Note: Modify E_List.cc
		int	E_List::operator!() 
		{ 
			return(head->getNext() != NULL); //modify from "== NULL" to "!=NULL", by MongTingWu
		}
Execute: $./LWHP.out
Input: Simulation times and E[td]
Output: Rh[Td] (Both mathematics and simulation analysis)
