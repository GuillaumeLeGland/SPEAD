function [success, OSMessage]=mycopyfile(src,dest)
    
%***********************************************************************************
%===================================================================================
% Use: [success,OSmessage]=mycopyfile(Source,Destination)
%===================================================================================
%  COPYFILE   Copy file or directory.
%     [SUCCESS,MESSAGE,MESSAGEID] = COPYFILE(SOURCE,DESTINATION,MODE) copies the
%     file or directory SOURCE to the new file or directory DESTINATION. Both
%     SOURCE and DESTINATION may be either an absolute pathname or a pathname
%     relative to the current directory. When the MODE is set, COPYFILE copies
%     SOURCE to DESTINATION, even when DESTINATION is read-only. The DESTINATION's
%     writable attribute state is preserved. See NOTE 1.
% 
%     [SUCCESS,MESSAGE,MESSAGEID] = COPYFILE(SOURCE) attempts to copy SOURCE to
%     the current directory.
% 
%     [SUCCESS,MESSAGE,MESSAGEID] = COPYFILE(SOURCE, DESTINATION) attempts to copy
%     SOURCE to DESTINATION. If SOURCE constitutes a directory or multiple files
%     and DESTINATION does not exist, COPYFILE attempts to create DESTINATION as a
%     directory and copy SOURCE to DESTINATION. If SOURCE constitutes a directory
%     or multiple files and DESTINATION exists as a directory, COPYFILE attempts
%     to copy SOURCE to DESTINATION. If SOURCE constitutes a directory or multiple
%     files and none of the above cases on DESTINATION applies, COPYFILE fails.
% 
%     [SUCCESS,MESSAGE,MESSAGEID] = COPYFILE(SOURCE,DESTINATION,'f') attempts to
%     copy SOURCE to DESTINATION, as above, even if DESTINATION is read-only. The
%     status of the writable attribute of DESTINATION will be preserved.
% 
%     INPUT PARAMETERS:
%         SOURCE:      1 x n string, defining the source file or directory.
%                      See NOTE 2 and 3.
%         DESTINATION: 1 x n string, defining destination file or directory.
%                      The default is the current directory. See NOTE 3.
%         MODE:        character scalar defining copy mode.
%                      'f' : force SOURCE to be written to DESTINATION. If omitted,
%                      COPYFILE respects the current writable status of DESTINATION.
%                      See NOTE 4.
% 
%     RETURN PARAMETERS:
%         SUCCESS:     logical scalar, defining the outcome of COPYFILE.
%                      1 : COPYFILE executed successfully.
%                      0 : an error occurred.
%         MESSAGE:     string, defining the error or warning message.
%                      empty string : COPYFILE executed successfully.
%                      message : an error or warning message, as applicable.
%         MESSAGEID:   string, defining the error or warning identifier.
%                      empty string : COPYFILE executed successfully.
%                      message id: the MATLAB error or warning message identifier
%                      (see ERROR, LASTERR, WARNING, LASTWARN).
% 
%     NOTE 1: Currently, SOURCE attributes are not preserved under copying on a
%             Windows platform. Except where otherwise stated, the rules of the
%             underlying system on the preservation of attributes are followed
%             when copying files and directories.
%     NOTE 2: The * wildcard, as a suffix to the last name or the extension to the
%             last name in a path string, is supported. Current behaviour of
%             COPYFILE differs between UNIX and Windows when using the wildcard *
%             or copying directories. See DOC COPYFILE for details.
%     NOTE 3: UNC paths are supported.
%     NOTE 4: 'writable' is being deprecated, but is still supported for backwards
%             compatibility.
% 
%     See also CD, DELETE, DIR, FILEATTRIB, MKDIR, MOVEFILE, RMDIR.
%===================================================================================
%     <http://neutron.risoe.dk/documentation/problems/>
%     <http://www.mathworks.com/support/solutions/en/data/1-1B5JY/index.html?solution=1-1B5JY>
%===================================================================================

if ispc
    [Status, OSMessage] = dos(['copy ' src ' '  dest]);
elseif isunix
    [Status, OSMessage] = unix(['cp -r ' src ' ' dest]);
end
success = ~Status;
