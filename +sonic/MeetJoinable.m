classdef (Abstract) MeetJoinable

    % Abstract parent class for all objects that can perform meet and join
    % operations. Just used to route these meets/joins to the appropriate
    % handling methods. 

    methods (Static, Hidden)

        function res = meetJoinWrapper(oper, handles, type_pairs, varargin)
        %% res = meetJoinWrapper(oper, handles, type_pairs, varargin)
        %
        %   Given a set of possible operations corresponding to certain
        %   datatypes, directs an arbitrary number of geometric objects to
        %   the appropriate operation handlers.
        %
        %   Inputs:
        %       - oper (1x1 string or char): String specifying the type of
        %         operation ('meet' or 'join') being conducted. Used for
        %         pretty-printing.
        %       - handles (nx1 cell array): Handler to call depending on
        %         object type. The i-th entry in `handles` corresponds to the
        %         i-th entry in `type_pairs`. 
        %       - type_pairs (nx1 cell array of 2x1 cell arrays): The list
        %         of datatypes that should trigger calling a certain handler.
        %         The i-th entry in `type_pairs` corresponds to the i-th 
        %         entry in `handles`. 
        %       - varargin (m sonic.MeetJoinable): A number of objects to
        %         either meet or join. 
        %
        %   Outputs:
        %       - res (1x1 sonic.MeetJoinable): The result of the chain of
        %         meet or join operations.
        %
        %   Last revised: 3/28/24
        %   Last author: Michael Krause

            % Check to ensure we have >= 2 inputs to meet/join:
            if length(varargin) < 2
                error('sonic:MeetJoinable:meetJoinWrapper:notEnoughObjects', ...
                      ['The %s operator requires at least ' ...
                       'two objects. Specify at least two objects to ' ...
                       'calculate the %s.'], oper, oper);
            end

            % Grab pairs and process:
            last_obj = varargin{1};
            for idx = 2:length(varargin)
               
                % Get the types of objects we're operating on:
                types = sonic.MeetJoinable.getTypes(last_obj, varargin{idx});

                % Check if the types match any of our type pairs:
                handle_idx = -1;
                flip_objs = false;
                for jdx = 1:length(type_pairs)
                    [do_match, do_flip] = ...
                        sonic.MeetJoinable.strPairMatchUnordered(...
                            types, type_pairs{jdx} ...
                        );
    
                    if do_match
                        handle_idx = jdx;
                        flip_objs = do_flip;
                        break;
                    end
    
                end

                % If we find a handle that matches, execute it:
                if handle_idx > 0
                    handle = handles{handle_idx};
                    if flip_objs
                        last_obj = handle(varargin{idx}, last_obj);
                    else
                        last_obj = handle(last_obj, varargin{idx});
                    end
                else

                    % For the user's convenience, we should pretty-print
                    % all possible options for them on error. This will
                    % require a bit of string manipulation to make it look
                    % nice:

                    oper_cases = '';
                    for jdx = 1:length(type_pairs)
                        oper_cases = [oper_cases ...
                            sprintf('\t - %s and %s\n', type_pairs{jdx}{1}, type_pairs{jdx}{2})];                        
                    end

                    % Fire off the actual error:
                    error('sonic:MeetJoinable:meetJoinWrapper:noHandlerFound', ...
                        ['Could not find a valid %s handler for the type ' ...
                        'combination of %s and %s. Valid combinations ' ...
                        'are:\n %s'], ...
                        oper, types{1}, types{2}, oper_cases);

                end

            end

            % And assign the output as the most recent output of one of our
            % meet/join handlers:
            res = last_obj;
        end

        function [do_match, do_flip] = strPairMatchUnordered(A, B)
        %% [do_match, do_flip] = strPairMatchUnordered(A, B)
        %
        %   Given two cell arrays, each of length two, this method checks
        %   if the two cell arrays contain the same contents, regardless of
        %   the order of the contents. It also indicates if they are in the
        %   same order or opposite order, if they do indeed match. 
        %
        %   Inputs:
        %       - A (2x1 cell): Cell array containing two strings 
        %       - B (2x1 cell): Cell array containing two strings
        %
        %   Outputs:
        %       - do_match (1x1 logical): True if the contents of A and B
        %         are the same (in the set-wise sense).
        %       - do_flip (1x1 logical): If do_match is true, this will be
        %         true if the ordering of items in A and B are flipped. It
        %         will be false in all other cases. 
        %
        %   Last revised: 3/28/24
        %   Last author: Michael Krause

            do_match = false;
            do_flip = false;

            if strcmp(A{1}, B{1}) && strcmp(A{2}, B{2})
                do_match = true;
                do_flip = false;
            elseif strcmp(A{1}, B{2}) && strcmp(A{2}, B{1})
                do_match = true;
                do_flip = true;
            end

        end

        function types = getTypes(varargin)
        %% types = getTypes(varargin)
        %
        %   Gets the types of an arbitrary number of objects passed in to
        %   this method. Given the application of this method, it also
        %   checks if each object is a subclass of sonic.MeetJoinable, and
        %   throws an error if not.
        %
        %   Inputs:
        %       - varargin (n): a variable number of objects
        %
        %   Outputs:
        %       - types (nx1 cell): the type of each object. Given
        %         conditions in the method, this will only return if it is a
        %         subclass of sonic.MeetJoinable, so the `sonic.` in the type
        %         designator is omitted in this vector. 
        %
        %   Last revised: 3/28/24
        %   Last author: Michael Krause

            types = cell(1, nargin);

            for idx = 1:length(types)

                % Check that this object is Meetable/Joinable.
                if ~isa(varargin{idx}, 'sonic.MeetJoinable')
                    error('sonic:MeetJoinable:getTypes:notMeetJoinable', ...
                        ['Object passed in to meet/join operation is ' ...
                        'not a subclass of sonic.MeetJoinable, and ' ...
                        'thus cannot be used for this operation.']);
                end

                % We now know we have a SONIC class, so once we get it,
                % let's chop off the 'sonic.' prefix for easier matching.
                raw_type = class(varargin{idx});
                split_type = split(raw_type, '.');
                types{idx} = split_type{2};
            end
            
        end

        function [obj_mult, obj_single] = verifyOneObjWithMultipleN(obj1, obj2)
        %% [obj_mult, obj_single] = verifyOneObjWithMultipleN(obj1, obj2)
        %
        %   Given two geometry objects which could contain multiple
        %   geometrical objects (such as Points or Lines objects),
        %   this method tests how many geometrical objects are contained
        %   within each object. If both `obj1` and `obj2` have n > 1, an
        %   error will be thrown. If both `obj1` and `obj2` have n = 1,
        %   then the `has_mult` flag will be set to false, indicating that
        %   neither object contains multiple geometrical instances.
        %   However, if exactly one of `obj1` or `obj2` has n > 1, while
        %   the other has n = 1, then the `has_mult` flag is set to true,
        %   and the corresponding object with multiple instances to
        %   assigned to `obj_mult` (likewise, the other is set to
        %   `obj_single`). This workflow enables quick downstream
        %   processing by meet/join methods that support many-to-one
        %   meets/joins.
        %
        %   Inputs:
        %       - obj1 (1x1 sonic.MeetJoinable): An object undergoing a
        %         meet/join that may or may not contain more than one
        %         geometrical instance.
        %       - obj2 (1x1 sonic.MeetJoinable): An object undergoing a
        %         meet/join that may or may not contain more than one
        %         geometrical instance.
        %
        %   Outputs:
        %       - obj_mult (1x1 sonic.MeetJoinable): If `has_mult` is true,
        %         this is the object (either `obj1` or `obj2`) with `n` > 1.
        %         If `has_mult` is false, this is [].
        %       - obj_single (1x1 sonic.MeetJoinable): If `has_mult` is 
        %         true, this is the object (either `obj1` or `obj2`) with 
        %         `n` = 1. If `has_mult` is false, this is [].
        %
        %   Last revised: 4/18/24
        %   Last author: Michael Krause

            n_gt1 = [obj1.n, obj2.n] > 1;

            if sum(n_gt1) > 1

                type1 = class(obj1);
                type2 = class(obj2);

                error('sonic:MeetJoinable:verifyOneObjWithMultipleN:tooManyObj', ...
                    ['When meeting or joining objects that can ' ...
                    'contain multiple instances of geometry (such as ' ...
                    'a Lines2/3 or Points2/3 object), only one object ' ...
                    'may contain multiple instances of geometry. ' ...
                    'Object of type %s contains %d instances, ' ...
                    'and object of type %s contains %d instances.'], ...
                    type1, obj1.n, type2, obj2.n);

            end

            % Check if we have any objects with multiple entries. If we do,
            % then return these objects accordingly. This allows for quick
            % processing by downstream methods. 
            if sum(n_gt1) > 0
                if n_gt1(1)
                    obj_mult = obj1;
                    obj_single = obj2;
                else
                    obj_mult = obj2;
                    obj_single = obj1;
                end
            else
                obj_mult = obj1; 
                obj_single = obj2; 
            end

        end
    
    end

end